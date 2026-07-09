#!/usr/bin/env python
"""Run an installed Python workflow script with concise user-facing errors."""
from __future__ import annotations

import os
import runpy
import sys
import traceback
from pathlib import Path


def prepend_active_env_bin_to_path() -> None:
    """Make subprocess calls use tools from the same conda environment."""
    env_bin = Path(sys.executable).resolve().parent
    old_path = os.environ.get("PATH", "")
    parts = old_path.split(os.pathsep) if old_path else []
    env_bin_str = str(env_bin)
    if not parts or parts[0] != env_bin_str:
        os.environ["PATH"] = os.pathsep.join([env_bin_str, *[p for p in parts if p != env_bin_str]])
    os.environ.setdefault("CONDA_PREFIX", str(env_bin.parent))


def main() -> int:
    if len(sys.argv) < 2:
        print("ERROR: missing script path.", file=sys.stderr)
        return 2
    script = Path(sys.argv[1]).resolve()
    if not script.is_file():
        print(f"ERROR: script not found: {script}", file=sys.stderr)
        return 2
    prepend_active_env_bin_to_path()
    sys.argv = [str(script)] + sys.argv[2:]
    try:
        runpy.run_path(str(script), run_name="__main__")
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        print(code, file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("ERROR: interrupted.", file=sys.stderr)
        return 130
    except Exception as exc:
        if os.environ.get("SPGRN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}:
            traceback.print_exc()
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
