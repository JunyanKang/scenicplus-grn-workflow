"""AutoZyme activation helpers for the SCENIC+ GRN workflow.

The helpers intentionally never install dependencies and never fail an analysis
because an AutoZyme patch family is unavailable. Installation is handled by the
workflow with no-deps mode; activation is auditable and can be disabled with
AUTOZYME_DISABLED=1.
"""
from __future__ import annotations

import os
import sys
from typing import Iterable


def activate_autozyme_python(families: Iterable[str] = ("scanpy",)) -> dict[str, str]:
    status: dict[str, str] = {}
    if os.environ.get("AUTOZYME_DISABLED", "0") == "1":
        status["autozyme"] = "disabled_by_env"
        print("[autozyme] disabled by AUTOZYME_DISABLED=1", file=sys.stderr)
        return status
    try:
        import autozyme  # type: ignore
    except Exception as exc:  # pragma: no cover - diagnostic path
        status["autozyme"] = f"not_available: {exc}"
        print(f"[autozyme] not available: {exc}", file=sys.stderr)
        return status

    for family in families:
        try:
            autozyme.activate(family)
            status[family] = "active"
        except Exception as exc:  # pragma: no cover - depends on upstream package
            status[family] = f"not_active: {exc}"
            print(f"[autozyme] could not activate {family}: {exc}", file=sys.stderr)
    try:
        full_status = autozyme.status()
        print(f"[autozyme] status: {full_status}", file=sys.stderr)
    except Exception:
        pass
    return status
