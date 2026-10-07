# -*- coding: utf-8 -*-
"""Part 2 후보 측정 -> CSV (보고서 ④ 달성표).

    python -m part2.search2               # 모든 후보
    python -m part2.search2 --only A
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from part2 import spec2
from part2.cell import CELLS

OUT = Path(__file__).resolve().parent.parent / "results"

COLS = [("c_store_fF", "C_STORE_fF", "-"),
        ("v_high_V", "V_high_V", "-"),
        ("dvbl_high_mV", "dVBL_high_mV", ">= 40"),
        ("dvbl_low_mV", "dVBL_low_mV", ">= 40"),
        ("retention_low_conservative_ms", "retention_low_ms", ">= 64"),
        ("retention_high_conservative_ms", "retention_high_ms", ">= 64"),
        ("retention_low_ms", "retention_low_trajectory_ms", "-"),
        ("retention_high_ms", "retention_high_trajectory_ms", "-"),
        ("eox_MV_per_cm", "Eox_MV_per_cm", "<= 5"),
        ("ediel_MV_per_cm", "Ediel_MV_per_cm", "<= 4")]
DEV = [("gate_length_um", "Lg_um"), ("oxide_thickness_nm", "tox_nm"),
       ("junction_depth_um", "xj_um"), ("silicon_thickness_um", "tSi_um"),
       ("body_doping_cm3", "NA_cm3"), ("sd_doping_cm3", "ND_cm3"),
       ("gate_material", "gate_metal"), ("cap_height_um", "cap_h_um"),
       ("cap_dielectric_nm", "cap_t_nm"), ("cap_pillar_nm", "pillar_nm"),
       ("cap_material", "cap_material")]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    keys = args.only or list(CELLS)

    rows = []
    for k in keys:
        r = spec2.measure(CELLS[k])
        rows.append((k, r))
        print(f"\n===== {k} =====\n{spec2.report(r)}", flush=True)

    OUT.mkdir(exist_ok=True)
    path = OUT / "part2_specs.csv"
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["cell"] + [d[1] for d in DEV] + [c[1] for c in COLS]
                   + [c[1] + "_pass" for c in COLS] + ["all_pass", "ileak_low_pA_per_um",
                                                       "ileak_high_pA_per_um", "range_violations", "seconds"])
        w.writerow(["기준"] + [""] * len(DEV) + [c[2] for c in COLS] + [""] * (len(COLS) + 5))
        for k, r in rows:
            d = r["device"]
            w.writerow([k] + [d.get(f, "") for f, _ in DEV]
                       + [r.get(c, "") for c, _, _ in COLS]
                       + ["PASS" if r["pass"].get(c) else "FAIL" for c, _, _ in COLS]
                       + [r.get("all_pass", ""), r.get("ileak_low_pA_per_um", ""),
                          r.get("ileak_high_pA_per_um", ""),
                          "; ".join(r.get("range_violations", [])),
                          round(float(r.get("seconds", 0) or 0), 1)])
    (OUT / "part2_specs.json").write_text(
        json.dumps({k: r for k, r in rows}, indent=2, default=str), encoding="utf-8")
    print(f"\nsaved: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
