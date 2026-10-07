# -*- coding: utf-8 -*-
"""Part 1 설계 변수. step2 의 Device 에 게이트 재료 두 항목만 더한다."""

from __future__ import annotations

from dataclasses import dataclass, asdict

from mosfet_tool.config import Device

from .gate_physics import GATE_METALS, work_function


@dataclass
class Part1Device(Device):
    """step2 Device + 게이트 재료.

    gate_material 이름이 유일한 진실 원천이다. 일함수 숫자를 따로 들고 있으면
    .devsim 파일에 들어가는 이름과 접점 방정식에 들어가는 숫자가 어긋날 수 있고,
    그건 그대로 감점으로 이어진다.
    """

    gate_material: str = "TaN"
    gate_metal_thickness_nm: float = 50.0   # 라벨 영역 두께. 물리적 효과 없음

    def __post_init__(self) -> None:
        work_function(self.gate_material)   # 표에 없는 이름이면 여기서 바로 막는다

    @property
    def gate_work_function_ev(self) -> float:
        return work_function(self.gate_material)

    def summary(self) -> dict:
        d = asdict(self)
        d["gate_work_function_ev"] = self.gate_work_function_ev
        d["Eox_MV_per_cm"] = 2.0 / (self.oxide_thickness_nm * 1e-7) / 1e6
        return d
