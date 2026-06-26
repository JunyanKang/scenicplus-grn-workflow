#!/usr/bin/env python
"""Run an installed Python workflow script with concise user-facing errors."""
from __future__ import annotations

import os
import runpy
import sys
import traceback
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 2:
        print("ERROR: missing script path.", file=sys.stderr)
        return 2
    script = Path(sys.argv[1]).resolve()
    if not script.is_file():
        print(f"ERROR: script not found: {script}", file=sys.stderr)
        return 2
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
