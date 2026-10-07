# -*- coding: utf-8 -*-
"""Part 1 후보 설계 모음. 탐색 기록표와 최종 제출 설계가 여기 한 곳에 모인다."""

from __future__ import annotations

from .device import Part1Device

# 배포 기본 소자 (비교 기준점)
BASELINE = Part1Device(
    gate_length_um=1.0, source_length_um=0.5, drain_length_um=0.5,
    oxide_thickness_nm=10.0, silicon_thickness_um=0.5, junction_depth_um=0.1,
    body_doping_cm3=1.0e16, sd_doping_cm3=1.0e19,
    gate_material="W",
)

# 후보 A — 짧은 채널 · 높은 기판 도핑 · TaN
CAND_A = Part1Device(
    gate_length_um=0.50, source_length_um=0.50, drain_length_um=0.50,
    oxide_thickness_nm=5.5, silicon_thickness_um=0.50, junction_depth_um=0.05,
    body_doping_cm3=7.5e16, sd_doping_cm3=1.0e20,
    gate_material="TaN",
)

# 후보 B — 긴 채널 · 얇은 산화막 · 두꺼운 실리콘 · W
CAND_B = Part1Device(
    gate_length_um=0.65, source_length_um=0.50, drain_length_um=0.50,
    oxide_thickness_nm=4.75, silicon_thickness_um=1.00, junction_depth_um=0.05,
    body_doping_cm3=1.8e16, sd_doping_cm3=1.0e20,
    gate_material="W",
)

# 최종 설계 — 후보 A 에서 산화막만 5.5 -> 5.0 nm 로 줄인 것.
# 탐색 10건 중 "조교가 다른 물리로 다시 돌렸을 때 흔들릴 수 있는 항목"의 여유가
# 가장 넓다. Vth 0.4505 V 는 허용 창(0.45 +- 0.05)의 정중앙이고, Ion 은 42% 여유,
# 고온 Ioff 는 2.3배 여유다. 가장 빡빡한 항목은 SS(예산의 93%)인데, SS 는
# C_dm/C_ox 로 정해지는 순수 정전 문제라 물리 모델이 달라져도 거의 움직이지 않는다.
FINAL_DEVICE = Part1Device(
    gate_length_um=0.50, source_length_um=0.50, drain_length_um=0.50,
    oxide_thickness_nm=5.0, silicon_thickness_um=0.50, junction_depth_um=0.05,
    body_doping_cm3=7.5e16, sd_doping_cm3=1.0e20,
    gate_material="TaN",
)

DESIGNS = {
    "baseline": BASELINE,
    "A": CAND_A,
    "B": CAND_B,
    "final": FINAL_DEVICE,
}

FINAL = "final"      # 제출 설계
