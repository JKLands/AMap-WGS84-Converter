# AMap GCJ-02 to WGS-84 GUI

高德 GCJ-02 点数据转换为 WGS-84 的 Windows 图形界面工具。

## 目录

- `amap_gdb_to_wgs84_gui.py`：核心源码
- `amap_gdb_to_wgs84_gui_v21.py`：v2.1 启动脚本，保留 FileGDB Integer64 字段
- `使用说明.txt`：使用说明

## 本地运行

需要 Python 3.10 或更高版本，并安装：

```powershell
python -m pip install -r requirements.txt
python amap_gdb_to_wgs84_gui_v21.py
```

## Windows 免安装版

请下载发布包中的 `网盘运行版`，保持 `AMap_GCJ02_to_WGS84_GUI.exe` 与 `_internal` 文件夹的相对位置，然后双击 EXE。

## 隐私提示

程序不会预设个人路径或账号信息。转换完成后生成的核验报告会记录用户选择的源路径和输出路径，请在公开分享前检查报告内容。
