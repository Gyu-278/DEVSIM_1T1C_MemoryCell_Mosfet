# -*- coding: utf-8 -*-
"""Ablation (보고서 ③). 최종 설계에서 파라미터를 하나씩만 바꿔 어느 지표가 움직이는지 본다.

    python -m part1.ablation --design A
"""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

from . import csvout, spec
from .designs import DESIGNS, FINAL

OUT = Path(__file__).resolve().parent.parent / "results"

# (필드, 표시 이름, 시도할 값들)
KNOBS = [
    ("oxide_thickness_nm", "tox [nm]", [4.0, 5.5, 8.0]),
    ("body_doping_cm3", "NA [cm^-3]", [1.0e16, 7.5e16, 1.0e17]),
    ("gate_length_um", "Lg [um]", [0.3, 0.5, 0.8]),
    ("silicon_thickness_um", "t_Si [um]", [0.3, 0.5, 1.0]),
    ("junction_depth_um", "xj [um]", [0.03, 0.05, 0.15]),
    ("sd_doping_cm3", "ND [cm^-3]", [1.0e19, 1.0e20, 1.0e21]),
]

COLS = ["vth", "ion_uA_per_um", "ioff_pA_per_um", "ss_mV_per_dec",
        "ioff_hot_pA_per_um", "dibl_mV_per_V", "body_effect_mV", "eox_MV_per_cm"]
HEAD = ["Vth", "Ion", "Ioff300", "SS", "Ioff398", "DIBL", "Body", "Eox"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default=FINAL, choices=sorted(DESIGNS))
    args = ap.parse_args()
    base = DESIGNS[args.design]

    rows = []
    print(f"reference = design {args.design}")
    ref = spec.measure(base)
    rows.append(("reference", ref))
    print(spec.report(ref), flush=True)

    for field, label, values in KNOBS:
        for v in values:
            if getattr(base, field) == v:
                continue
            dev = dataclasses.replace(base, **{field: v})
            res = spec.measure(dev)
            tag = f"{label}={v:g}"
            rows.append((tag, res))
            deltas = []
            for c, h in zip(COLS, HEAD):
                a, b = ref.get(c), res.get(c)
                if isinstance(a, float) and isinstance(b, float) and a == a and b == b and a != 0:
                    deltas.append(f"{h} {(b/a-1)*100:+.0f}%")
            print(f"{tag:22s} " + "  ".join(deltas) +
                  ("   ALL PASS" if res["all_pass"] else "   fail: " +
                   ",".join(k for k, ok in res["pass"].items() if not ok)), flush=True)

    OUT.mkdir(exist_ok=True)
    (OUT / f"ablation_{args.design}.json").write_text(
        json.dumps({n: r for n, r in rows}, indent=2, default=str), encoding="utf-8")
    csv_path = csvout.write_specs(rows, OUT / f"ablation_{args.design}.csv")
    print(f"\nsaved: {OUT / f'ablation_{args.design}.json'} , {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
