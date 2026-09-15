#!/usr/bin/env python3
"""v2.1 launcher with FileGDB Integer64 preservation enabled."""

from __future__ import annotations

import amap_gdb_to_wgs84_gui as core


_pyogrio_write = core.write


def write_with_filegdb_compatibility(*args, **kwargs):
    """Keep Integer64 fields when writing modern OpenFileGDB outputs."""
    if kwargs.get("driver") == "OpenFileGDB":
        options = dict(kwargs.get("layer_options") or {})
        options.setdefault("TARGET_ARCGIS_VERSION", "ARCGIS_PRO_3_2_OR_LATER")
        kwargs["layer_options"] = options
    return _pyogrio_write(*args, **kwargs)


core.write = write_with_filegdb_compatibility
core.APP_VERSION = "2.1"


if __name__ == "__main__":
    raise SystemExit(core.launch_gui())
