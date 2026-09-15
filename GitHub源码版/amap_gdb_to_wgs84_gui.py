#!/usr/bin/env python3
"""通用高德 GCJ-02 点数据转 WGS-84 图形界面程序。

读取用户选择的 Esri FileGDB，转换其中全部 Point 图层；源数据不会被修改。
支持输出新的 FileGDB（.gdb）或 GeoPackage（.gpkg）。
"""

from __future__ import annotations

import datetime as dt
import json
import math
import os
import shutil
import struct
import sys
import threading
import traceback
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pyogrio
from pyogrio.raw import read, write


APP_NAME = "高德 GCJ-02 → WGS-84 通用转换工具"
APP_VERSION = "2.0"
PI = math.pi
SEMI_MAJOR_AXIS = 6378245.0
ECCENTRICITY_SQUARED = 0.00669342162296594323
ProgressCallback = Callable[[str], None]


def in_gcj_region(lng: float, lat: float) -> bool:
    """Conservative bounding box used to avoid modifying obvious overseas data."""
    return 72.004 <= lng <= 137.8347 and 0.8293 <= lat <= 55.8271


def transform_latitude(x: float, y: float) -> float:
    value = (
        -100.0
        + 2.0 * x
        + 3.0 * y
        + 0.2 * y * y
        + 0.1 * x * y
        + 0.2 * math.sqrt(abs(x))
    )
    value += (
        20.0 * math.sin(6.0 * x * PI) + 20.0 * math.sin(2.0 * x * PI)
    ) * 2.0 / 3.0
    value += (
        20.0 * math.sin(y * PI) + 40.0 * math.sin(y / 3.0 * PI)
    ) * 2.0 / 3.0
    value += (
        160.0 * math.sin(y / 12.0 * PI) + 320.0 * math.sin(y * PI / 30.0)
    ) * 2.0 / 3.0
    return value


def transform_longitude(x: float, y: float) -> float:
    value = (
        300.0
        + x
        + 2.0 * y
        + 0.1 * x * x
        + 0.1 * x * y
        + 0.1 * math.sqrt(abs(x))
    )
    value += (
        20.0 * math.sin(6.0 * x * PI) + 20.0 * math.sin(2.0 * x * PI)
    ) * 2.0 / 3.0
    value += (
        20.0 * math.sin(x * PI) + 40.0 * math.sin(x / 3.0 * PI)
    ) * 2.0 / 3.0
    value += (
        150.0 * math.sin(x / 12.0 * PI) + 300.0 * math.sin(x / 30.0 * PI)
    ) * 2.0 / 3.0
    return value


def wgs84_to_gcj02(lng: float, lat: float) -> tuple[float, float]:
    delta_lat = transform_latitude(lng - 105.0, lat - 35.0)
    delta_lng = transform_longitude(lng - 105.0, lat - 35.0)
    rad_lat = lat / 180.0 * PI
    magic = math.sin(rad_lat)
    magic = 1.0 - ECCENTRICITY_SQUARED * magic * magic
    sqrt_magic = math.sqrt(magic)
    delta_lat = (
        delta_lat
        * 180.0
        / (
            (SEMI_MAJOR_AXIS * (1.0 - ECCENTRICITY_SQUARED))
            / (magic * sqrt_magic)
            * PI
        )
    )
    delta_lng = (
        delta_lng
        * 180.0
        / (SEMI_MAJOR_AXIS / sqrt_magic * math.cos(rad_lat) * PI)
    )
    return lng + delta_lng, lat + delta_lat


def gcj02_to_wgs84(
    gcj_lng: float,
    gcj_lat: float,
    iterations: int = 10,
    tolerance: float = 1e-9,
) -> tuple[float, float]:
    if not in_gcj_region(gcj_lng, gcj_lat):
        return float(gcj_lng), float(gcj_lat)
    wgs_lng = float(gcj_lng)
    wgs_lat = float(gcj_lat)
    for _ in range(iterations):
        test_lng, test_lat = wgs84_to_gcj02(wgs_lng, wgs_lat)
        delta_lng = test_lng - gcj_lng
        delta_lat = test_lat - gcj_lat
        wgs_lng -= delta_lng
        wgs_lat -= delta_lat
        if abs(delta_lng) <= tolerance and abs(delta_lat) <= tolerance:
            break
    return wgs_lng, wgs_lat


def point_layout(wkb: bytes) -> tuple[str, int]:
    if len(wkb) < 21:
        raise ValueError("点几何 WKB 长度不足")
    endian = "<" if wkb[0] == 1 else ">" if wkb[0] == 0 else None
    if endian is None:
        raise ValueError("无法识别 WKB 字节序")
    type_code = struct.unpack_from(endian + "I", wkb, 1)[0]
    base_type = (type_code & 0x1FFFFFFF) % 1000
    if base_type != 1:
        raise ValueError(f"仅支持 Point 几何，收到 WKB 类型 {type_code}")
    coordinate_offset = 9 if type_code & 0x20000000 else 5
    if len(wkb) < coordinate_offset + 16:
        raise ValueError("点几何 WKB 坐标数据不完整")
    return endian, coordinate_offset


def point_from_wkb(wkb: bytes | None) -> tuple[float, float] | None:
    if wkb is None:
        return None
    endian, offset = point_layout(wkb)
    lng, lat = struct.unpack_from(endian + "dd", wkb, offset)
    if not (math.isfinite(lng) and math.isfinite(lat)):
        return None
    return lng, lat


def convert_point_wkb(
    wkb: bytes | None,
) -> tuple[bytes | None, tuple[float, float] | None, bool]:
    if wkb is None:
        return None, None, False
    endian, offset = point_layout(wkb)
    gcj_lng, gcj_lat = struct.unpack_from(endian + "dd", wkb, offset)
    if not (math.isfinite(gcj_lng) and math.isfinite(gcj_lat)):
        return wkb, None, False
    is_converted = in_gcj_region(gcj_lng, gcj_lat)
    wgs_lng, wgs_lat = gcj02_to_wgs84(gcj_lng, gcj_lat)
    output = bytearray(wkb)
    struct.pack_into(endian + "dd", output, offset, wgs_lng, wgs_lat)
    return bytes(output), (wgs_lng, wgs_lat), is_converted


def as_float_array(values: Iterable[Any]) -> np.ndarray:
    result: list[float] = []
    for value in values:
        if value is None or value == "":
            result.append(np.nan)
            continue
        try:
            result.append(float(value))
        except (TypeError, ValueError):
            result.append(np.nan)
    return np.asarray(result, dtype="float64")


def unique_field_name(preferred: str, fields: list[str]) -> str:
    existing = {field.casefold() for field in fields}
    if preferred.casefold() not in existing:
        return preferred
    index = 2
    while f"{preferred}_{index}".casefold() in existing:
        index += 1
    return f"{preferred}_{index}"


def find_field(fields: list[str], candidates: tuple[str, ...]) -> int | None:
    lookup = {field.casefold(): index for index, field in enumerate(fields)}
    for candidate in candidates:
        if candidate.casefold() in lookup:
            return lookup[candidate.casefold()]
    return None


def transform_attribute_coordinates(
    fields: list[str], field_data: list[np.ndarray]
) -> tuple[list[str], list[np.ndarray], dict[str, Any]]:
    output_fields = list(fields)
    output_data = list(field_data)
    lon_index = find_field(fields, ("lon", "lng", "longitude", "经度"))
    lat_index = find_field(fields, ("lat", "latitude", "纬度"))
    stats: dict[str, Any] = {
        "coordinate_fields_found": lon_index is not None and lat_index is not None,
        "attribute_coordinates_transformed": 0,
    }
    if lon_index is None or lat_index is None:
        return output_fields, output_data, stats

    source_lng = as_float_array(output_data[lon_index])
    source_lat = as_float_array(output_data[lat_index])
    result_lng = source_lng.copy()
    result_lat = source_lat.copy()
    valid = np.isfinite(source_lng) & np.isfinite(source_lat)
    converted = 0
    outside = 0
    for row_index in np.flatnonzero(valid):
        if in_gcj_region(source_lng[row_index], source_lat[row_index]):
            converted += 1
        else:
            outside += 1
        result_lng[row_index], result_lat[row_index] = gcj02_to_wgs84(
            source_lng[row_index], source_lat[row_index]
        )

    source_lon_name = unique_field_name("gcj02_lon", output_fields)
    output_fields.append(source_lon_name)
    output_data.append(source_lng)
    source_lat_name = unique_field_name("gcj02_lat", output_fields)
    output_fields.append(source_lat_name)
    output_data.append(source_lat)
    output_data[lon_index] = result_lng
    output_data[lat_index] = result_lat
    stats.update(
        {
            "longitude_field": fields[lon_index],
            "latitude_field": fields[lat_index],
            "original_longitude_field": source_lon_name,
            "original_latitude_field": source_lat_name,
            "attribute_coordinates_transformed": converted,
            "attribute_coordinates_outside_gcj_region_unchanged": outside,
            "attribute_coordinates_missing": int((~valid).sum()),
        }
    )
    return output_fields, output_data, stats


def calculate_bounds(points: list[tuple[float, float]]) -> list[float] | None:
    if not points:
        return None
    return [
        min(point[0] for point in points),
        min(point[1] for point in points),
        max(point[0] for point in points),
        max(point[1] for point in points),
    ]


def convert_layer(
    source: Path,
    output: Path,
    layer_name: str,
    driver: str,
) -> dict[str, Any]:
    meta, _fids, geometry, field_data = read(
        source, layer=layer_name, return_fids=True, force_2d=False
    )
    geometry_type = str(meta.get("geometry_type") or "")
    if not geometry_type.startswith("Point"):
        raise ValueError(f"图层 {layer_name!r} 不是点图层：{geometry_type}")

    source_points: list[tuple[float, float]] = []
    output_points: list[tuple[float, float]] = []
    converted_geometry = np.empty(len(geometry), dtype=object)
    converted_count = 0
    outside_count = 0
    empty_count = 0
    for index, item in enumerate(geometry):
        source_point = point_from_wkb(item)
        converted_wkb, output_point, was_converted = convert_point_wkb(item)
        converted_geometry[index] = converted_wkb
        if source_point is not None:
            source_points.append(source_point)
        if output_point is not None:
            output_points.append(output_point)
        if output_point is None:
            empty_count += 1
        elif was_converted:
            converted_count += 1
        else:
            outside_count += 1

    fields = [str(field) for field in meta["fields"]]
    output_fields, output_data, attribute_stats = transform_attribute_coordinates(
        fields, list(field_data)
    )
    write(
        output,
        converted_geometry,
        output_data,
        output_fields,
        layer=layer_name,
        driver=driver,
        geometry_type=geometry_type,
        crs="EPSG:4326",
        encoding="UTF-8",
        append=False,
        layer_metadata={
            "coordinate_conversion": "GCJ-02 to WGS-84 iterative numerical inverse",
            "source_crs_assumption": "GCJ-02",
        },
    )

    samples: list[dict[str, list[float]]] = []
    for source_point, output_point in zip(source_points, output_points):
        samples.append(
            {
                "gcj02": [round(source_point[0], 8), round(source_point[1], 8)],
                "wgs84": [round(output_point[0], 8), round(output_point[1], 8)],
            }
        )
        if len(samples) == 5:
            break

    return {
        "layer": layer_name,
        "features": len(geometry),
        "source_geometry_type": geometry_type,
        "source_crs_metadata": meta.get("crs"),
        "output_crs": "EPSG:4326",
        "geometry_coordinates_transformed": converted_count,
        "geometry_coordinates_outside_gcj_region_unchanged": outside_count,
        "null_or_empty_geometries": empty_count,
        "source_bounds": calculate_bounds(source_points),
        "output_bounds": calculate_bounds(output_points),
        "sample_coordinates": samples,
        **attribute_stats,
    }


def output_driver(output: Path) -> str:
    suffix = output.suffix.casefold()
    if suffix == ".gdb":
        return "OpenFileGDB"
    if suffix == ".gpkg":
        return "GPKG"
    raise ValueError("输出路径必须以 .gdb 或 .gpkg 结尾")


def stamped_sibling(path: Path, label: str) -> Path:
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    return path.with_name(f"{path.stem}.{label}_{stamp}{path.suffix}")


def convert_dataset(
    source: Path,
    output: Path,
    overwrite: bool = False,
    progress: ProgressCallback | None = None,
) -> dict[str, Any]:
    source = source.expanduser().resolve()
    output = output.expanduser().resolve()
    notify = progress or (lambda _message: None)
    if not source.exists() or not source.is_dir() or source.suffix.casefold() != ".gdb":
        raise FileNotFoundError("请选择有效的源 FileGDB 文件夹（扩展名为 .gdb）")
    if source == output:
        raise ValueError("输出路径不能与源 GDB 相同")
    driver = output_driver(output)
    if output.exists() and not overwrite:
        raise FileExistsError("输出已存在。请选择新位置，或勾选“备份并替换已有输出”。")

    available = [(str(name), str(kind)) for name, kind in pyogrio.list_layers(source)]
    point_layers = [name for name, kind in available if kind.startswith("Point")]
    if not point_layers:
        raise ValueError("源 GDB 中没有 Point 点图层")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = stamped_sibling(output, "partial")
    if temporary.exists():
        raise FileExistsError(f"临时输出已存在：{temporary}")

    layer_reports: list[dict[str, Any]] = []
    try:
        for index, layer_name in enumerate(point_layers, start=1):
            notify(f"[{index}/{len(point_layers)}] 正在转换图层：{layer_name}")
            layer_reports.append(convert_layer(source, temporary, layer_name, driver))

        backup_output = None
        if output.exists():
            backup_output = stamped_sibling(output, "backup")
            shutil.move(str(output), str(backup_output))
        shutil.move(str(temporary), str(output))
    except Exception:
        if temporary.exists():
            failed = stamped_sibling(output, "failed")
            shutil.move(str(temporary), str(failed))
            notify(f"失败的临时输出已保留：{failed}")
        raise

    report_path = output.with_name(f"{output.stem}_conversion_report.json")
    backup_report = None
    if report_path.exists():
        backup_report = stamped_sibling(report_path, "backup")
        shutil.move(str(report_path), str(backup_report))
    report = {
        "application": APP_NAME,
        "version": APP_VERSION,
        "created_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "source": str(source),
        "output": str(output),
        "output_driver": driver,
        "method": "GCJ-02 to WGS-84 iterative numerical inverse",
        "official_amap_reverse_api": False,
        "source_dataset_unchanged": True,
        "backup_of_previous_output": str(backup_output) if backup_output else None,
        "backup_of_previous_report": str(backup_report) if backup_report else None,
        "layers": layer_reports,
    }
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    report["report_path"] = str(report_path)
    return report


def launch_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    try:
        import ctypes

        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass

    class ConverterApp:
        def __init__(self, root: tk.Tk) -> None:
            self.root = root
            self.root.title(f"{APP_NAME} v{APP_VERSION}")
            self.root.geometry("820x570")
            self.root.minsize(720, 520)
            self.source_var = tk.StringVar()
            self.output_var = tk.StringVar()
            self.overwrite_var = tk.BooleanVar(value=False)
            self.last_output: Path | None = None
            self.build_ui(ttk)

        def build_ui(self, ttk_module: Any) -> None:
            style = ttk_module.Style()
            if "vista" in style.theme_names():
                style.theme_use("vista")
            style.configure("Title.TLabel", font=("Microsoft YaHei UI", 16, "bold"))
            style.configure("Hint.TLabel", foreground="#555555")

            outer = ttk_module.Frame(self.root, padding=22)
            outer.pack(fill="both", expand=True)
            outer.columnconfigure(1, weight=1)
            outer.rowconfigure(7, weight=1)

            ttk_module.Label(outer, text=APP_NAME, style="Title.TLabel").grid(
                row=0, column=0, columnspan=3, sticky="w", pady=(0, 6)
            )
            ttk_module.Label(
                outer,
                text="适用于已确认是高德 GCJ-02 的 POI 点数据；源数据不会被修改。",
                style="Hint.TLabel",
            ).grid(row=1, column=0, columnspan=3, sticky="w", pady=(0, 20))

            ttk_module.Label(outer, text="源 FileGDB：").grid(
                row=2, column=0, sticky="w", padx=(0, 10)
            )
            self.source_entry = ttk_module.Entry(outer, textvariable=self.source_var)
            self.source_entry.grid(row=2, column=1, sticky="ew")
            self.source_button = ttk_module.Button(
                outer, text="选择源文件夹…", command=self.choose_source
            )
            self.source_button.grid(row=2, column=2, padx=(10, 0))

            ttk_module.Label(outer, text="输出位置：").grid(
                row=3, column=0, sticky="w", padx=(0, 10), pady=(14, 0)
            )
            self.output_entry = ttk_module.Entry(outer, textvariable=self.output_var)
            self.output_entry.grid(row=3, column=1, sticky="ew", pady=(14, 0))
            self.output_button = ttk_module.Button(
                outer, text="选择输出…", command=self.choose_output
            )
            self.output_button.grid(row=3, column=2, padx=(10, 0), pady=(14, 0))

            ttk_module.Checkbutton(
                outer,
                text="备份并替换已有输出",
                variable=self.overwrite_var,
            ).grid(row=4, column=1, sticky="w", pady=(12, 0))

            action_frame = ttk_module.Frame(outer)
            action_frame.grid(row=5, column=0, columnspan=3, sticky="ew", pady=(18, 12))
            self.convert_button = ttk_module.Button(
                action_frame, text="开始转换", command=self.start_conversion
            )
            self.convert_button.pack(side="left")
            self.open_button = ttk_module.Button(
                action_frame,
                text="打开输出文件夹",
                command=self.open_output_folder,
                state="disabled",
            )
            self.open_button.pack(side="left", padx=(10, 0))
            self.progressbar = ttk_module.Progressbar(action_frame, mode="indeterminate")
            self.progressbar.pack(side="right", fill="x", expand=True, padx=(24, 0))

            ttk_module.Label(outer, text="运行日志：").grid(
                row=6, column=0, columnspan=3, sticky="w", pady=(0, 6)
            )
            log_frame = ttk_module.Frame(outer)
            log_frame.grid(row=7, column=0, columnspan=3, sticky="nsew")
            log_frame.columnconfigure(0, weight=1)
            log_frame.rowconfigure(0, weight=1)
            self.log_text = tk.Text(
                log_frame,
                height=13,
                wrap="word",
                state="disabled",
                font=("Microsoft YaHei UI", 9),
            )
            scrollbar = ttk_module.Scrollbar(
                log_frame, orient="vertical", command=self.log_text.yview
            )
            self.log_text.configure(yscrollcommand=scrollbar.set)
            self.log_text.grid(row=0, column=0, sticky="nsew")
            scrollbar.grid(row=0, column=1, sticky="ns")
            self.log("请选择源 .gdb 文件夹和输出位置。")

        def choose_source(self) -> None:
            selected = filedialog.askdirectory(
                parent=self.root,
                title="选择源 FileGDB 文件夹（.gdb）",
                mustexist=True,
            )
            if not selected:
                return
            source = Path(selected)
            if source.suffix.casefold() != ".gdb":
                messagebox.showwarning(
                    "选择无效", "请选择扩展名为 .gdb 的 FileGDB 文件夹。", parent=self.root
                )
                return
            self.source_var.set(str(source))
            self.output_var.set(str(source.with_name(f"{source.stem}_wgs84.gdb")))
            self.log(f"已选择源文件：{source}")

        def choose_output(self) -> None:
            source_text = self.source_var.get().strip()
            source = Path(source_text) if source_text else Path.cwd() / "poi.gdb"
            selected = filedialog.asksaveasfilename(
                parent=self.root,
                title="选择输出位置",
                initialdir=str(source.parent),
                initialfile=f"{source.stem}_wgs84.gdb",
                defaultextension=".gdb",
                filetypes=[
                    ("Esri FileGDB", "*.gdb"),
                    ("GeoPackage", "*.gpkg"),
                    ("所有文件", "*.*"),
                ],
            )
            if selected:
                self.output_var.set(selected)
                self.log(f"输出将写入：{selected}")

        def log(self, message: str) -> None:
            timestamp = dt.datetime.now().strftime("%H:%M:%S")
            self.log_text.configure(state="normal")
            self.log_text.insert("end", f"[{timestamp}] {message}\n")
            self.log_text.see("end")
            self.log_text.configure(state="disabled")

        def log_from_worker(self, message: str) -> None:
            self.root.after(0, self.log, message)

        def set_busy(self, busy: bool) -> None:
            state = "disabled" if busy else "normal"
            for widget in (
                self.source_button,
                self.output_button,
                self.convert_button,
                self.source_entry,
                self.output_entry,
            ):
                widget.configure(state=state)
            if busy:
                self.progressbar.start(12)
            else:
                self.progressbar.stop()

        def start_conversion(self) -> None:
            source_text = self.source_var.get().strip().strip('"')
            output_text = self.output_var.get().strip().strip('"')
            if not source_text or not output_text:
                messagebox.showwarning(
                    "信息不完整", "请先选择源 FileGDB 和输出位置。", parent=self.root
                )
                return
            source = Path(source_text)
            output = Path(output_text)
            if output.suffix.casefold() not in {".gdb", ".gpkg"}:
                messagebox.showwarning(
                    "输出格式无效", "输出名称必须以 .gdb 或 .gpkg 结尾。", parent=self.root
                )
                return
            self.set_busy(True)
            self.open_button.configure(state="disabled")
            self.log("开始读取并转换点图层……")
            worker = threading.Thread(
                target=self.run_conversion,
                args=(source, output, self.overwrite_var.get()),
                daemon=True,
            )
            worker.start()

        def run_conversion(self, source: Path, output: Path, overwrite: bool) -> None:
            try:
                report = convert_dataset(
                    source, output, overwrite=overwrite, progress=self.log_from_worker
                )
            except Exception as exc:
                details = traceback.format_exc()
                self.root.after(0, self.conversion_failed, str(exc), details)
                return
            self.root.after(0, self.conversion_finished, report)

        def conversion_failed(self, message: str, details: str) -> None:
            self.set_busy(False)
            self.log(f"转换失败：{message}")
            if os.environ.get("AMAP_CONVERTER_DEBUG") == "1":
                self.log(details)
            messagebox.showerror("转换失败", message, parent=self.root)

        def conversion_finished(self, report: dict[str, Any]) -> None:
            self.set_busy(False)
            self.last_output = Path(report["output"])
            self.open_button.configure(state="normal")
            features = sum(layer["features"] for layer in report["layers"])
            self.log(f"转换完成：{features:,} 条记录")
            self.log(f"结果文件：{report['output']}")
            self.log(f"核验报告：{report['report_path']}")
            messagebox.showinfo(
                "转换完成",
                f"已转换 {features:,} 条记录。\n\n结果：{report['output']}",
                parent=self.root,
            )

        def open_output_folder(self) -> None:
            if self.last_output is not None:
                os.startfile(str(self.last_output.parent))

    root = tk.Tk()
    ConverterApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(launch_gui())
