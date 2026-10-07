# -*- coding: utf-8 -*-
"""과제가 정한 설계 파라미터 허용 범위 검사.

공식 검사기(check_structure_file.py)는 영역·접점·재료 이름만 본다. 도핑이나 치수가
허용 범위 안인지는 **보지 않는다**. 그런데 과제 27쪽은 조교가 제출 파일의 메시
좌표와 도핑 값에서 Lg·tox·NA 를 다시 계산해 범위를 확인한다고 적어 두었다.
즉 검사기는 통과하고 조교에게서 걸리는 구간이 존재한다. 그 구멍을 여기서 막는다.
"""

from __future__ import annotations

# (필드, 표시 이름, 하한, 상한, 단위). None 은 제한 없음.
RANGES = [
    ("gate_length_um", "Lg", 0.3, None, "um"),
    ("oxide_thickness_nm", "tox", 4.0, None, "nm"),
    ("silicon_thickness_um", "t_Si", 0.3, 2.0, "um"),
    ("junction_depth_um", "xj", 0.02, 0.25, "um"),
    ("body_doping_cm3", "NA", 1.0e14, 1.0e17, "cm^-3"),
    ("sd_doping_cm3", "ND", 1.0e18, 1.0e21, "cm^-3"),
    ("source_length_um", "S length", 0.2, 1.0, "um"),
    ("drain_length_um", "D length", 0.2, 1.0, "um"),
]

# Part 2 캡 (과제 34쪽)
CAP_RANGES = [
    ("cap_height_um", "cap h", 0.2, 1.5, "um"),
    ("cap_dielectric_nm", "cap t", 3.0, 10.0, "nm"),
    ("cap_pillar_nm", "pillar w", 20.0, 80.0, "nm"),
]


def check(device, v_wl: float = 2.5, cap: bool = True) -> list[str]:
    """위반 목록을 돌려준다. 빈 리스트면 통과."""
    bad = []
    rows = list(RANGES) + (CAP_RANGES if cap else [])
    for field, label, lo, hi, unit in rows:
        if not hasattr(device, field):
            continue
        v = getattr(device, field)
        if lo is not None and v < lo:
            bad.append(f"{label} = {v:g} {unit} < 하한 {lo:g}")
        if hi is not None and v > hi:
            bad.append(f"{label} = {v:g} {unit} > 상한 {hi:g}")

    # C_STORE 상한 (2026-09-27 개정 스펙). 해석식으로 미리 거른다.
    if cap and hasattr(device, "c_store_closed_form_fF"):
        c = device.c_store_closed_form_fF
        if c > 20.0:
            bad.append(f"C_STORE(해석식) = {c:.2f} fF > 상한 20")

    # 전계 제약
    e_ox = v_wl / (device.oxide_thickness_nm * 1e-7) / 1e6
    if e_ox > 5.0:
        bad.append(f"E_ox = {e_ox:.3f} MV/cm > 5 (V_WL {v_wl} V, tox {device.oxide_thickness_nm} nm)")
    if cap and hasattr(device, "cap_dielectric_nm"):
        e_d = 1.0 / (device.cap_dielectric_nm * 1e-7) / 1e6
        if e_d > 4.0:
            bad.append(f"유전막 전계 = {e_d:.3f} MV/cm > 4 (t {device.cap_dielectric_nm} nm)")
    return bad


def assert_ok(device, v_wl: float = 2.5, cap: bool = True) -> None:
    bad = check(device, v_wl, cap)
    if bad:
        raise ValueError("과제 허용 범위 위반:\n  - " + "\n  - ".join(bad))
