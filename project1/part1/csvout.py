# -*- coding: utf-8 -*-
"""측정 결과를 엑셀에서 바로 열리는 CSV 로 쓴다.

한글이 깨지지 않도록 utf-8-sig(BOM) 로 저장한다. 엑셀이 BOM 을 보고 UTF-8 로 읽는다.
"""

from __future__ import annotations

import csv
from pathlib import Path

# (결과 키, CSV 열 이름, 기준 문구)
SPEC_COLS = [
    ("vth", "Vth_V", "0.45 +- 0.05"),
    ("ion_uA_per_um", "Ion_uA_per_um", ">= 450"),
    ("ioff_pA_per_um", "Ioff300K_pA_per_um", "<= 1"),
    ("ss_mV_per_dec", "SS_mV_per_dec", "<= 75"),
    ("ioff_hot_pA_per_um", "Ioff398K_pA_per_um", "<= 100"),
    ("dibl_mV_per_V", "DIBL_mV_per_V", "<= 30"),
    ("body_effect_mV", "BodyEffect_mV", "<= 80"),
    ("eox_MV_per_cm", "Eox_MV_per_cm", "<= 5"),
]

DEV_COLS = [
    ("gate_length_um", "Lg_um"),
    ("oxide_thickness_nm", "tox_nm"),
    ("junction_depth_um", "xj_um"),
    ("silicon_thickness_um", "tSi_um"),
    ("source_length_um", "Ls_um"),
    ("drain_length_um", "Ld_um"),
    ("body_doping_cm3", "NA_cm3"),
    ("sd_doping_cm3", "ND_cm3"),
    ("gate_material", "gate_metal"),
    ("gate_work_function_ev", "gate_EWF_eV"),
    ("temperature_k", "T_K"),
]


def write_specs(rows, path) -> Path:
    """[(이름, 결과dict), ...] -> CSV 한 장. 설계 변수 + 8항목 + 합격 여부."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (["variant"] + [d[1] for d in DEV_COLS] + [c[1] for c in SPEC_COLS]
              + [c[1] + "_pass" for c in SPEC_COLS]
              + ["all_pass", "trustworthy", "vth_bracket_V", "ss_points",
                 "ramp_retries", "seconds", "errors"])
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        # 기준값을 두 번째 줄에 같이 적어 둔다 (보고서에 그대로 붙일 수 있게)
        w.writerow(["기준"] + [""] * len(DEV_COLS) + [c[2] for c in SPEC_COLS]
                   + [""] * len(SPEC_COLS) + [""] * 7)
        for name, r in rows:
            dev = r.get("device", {})
            line = [name]
            line += [dev.get(k, "") for k, _ in DEV_COLS]
            line += [r.get(k, "") for k, _, _ in SPEC_COLS]
            line += ["PASS" if r.get("pass", {}).get(k) else "FAIL" for k, _, _ in SPEC_COLS]
            line += [r.get("all_pass", ""), r.get("trustworthy", ""),
                     r.get("vth_bracket_V", ""), r.get("ss_points", ""),
                     r.get("retries", ""), round(float(r.get("seconds", 0) or 0), 1),
                     str(r.get("errors", "")) or ""]
            w.writerow(line)
    return path


def write_curve(vgs, ids, path, label="") -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Vg_V", "Id_A_per_um", "bias"])
        for v, i in sorted(zip(vgs, ids)):
            w.writerow([v, i, label])
    return path
