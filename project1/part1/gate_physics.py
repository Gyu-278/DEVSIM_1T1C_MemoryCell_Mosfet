# -*- coding: utf-8 -*-
"""[추가 물리 2] 게이트 금속 일함수.

배포 코드의 CreateOxideContact 는 게이트 노드를 그냥 `Potential - V_gate = 0` 으로
묶는다. simple_physics.py 76 행에 "# in the future, worry about workfunction" 이라는
주석이 그 함수 바로 위에 달려 있다. 즉 일함수 항이 아예 없다.

이 코드의 전위 기준은 진성 실리콘이다 (IntrinsicElectrons = n_i*exp(Potential/V_t)).
그래서 배포 코드는 일함수가 phi_i(= 4.61 eV) 인 게이트, 곧 midgap 게이트를 쓴 셈이다.
과제 28쪽의 금속 표를 쓰려면 기준 phi_i 로부터의 차이만큼 게이트 노드 전위를
내려 주어야 한다:

    Potential = V_gate - (phi_m - phi_i)

일함수가 크면 게이트 노드 전위가 내려가고 Vth 가 올라간다. 이 관계는 근사가 아니라
게이트 전압축을 평행이동시키는 항등식이어서 dVth/dphi_m = 1.000 V/eV 로 정확하다.
따라서 금속을 바꿀 때 SS · DIBL · body effect 는 전혀 변하지 않고,
고정된 절대 Vg 에서 읽는 Ion · Ioff · 고온 Ioff 만 변한다.
"""

from __future__ import annotations

from devsim import contact_equation
from devsim.python_packages.model_create import (
    CreateContactNodeModel,
    CreateContactNodeModelDerivative,
    CreateEdgeModel,
    CreateEdgeModelDerivatives,
    InEdgeModelList,
)
from devsim.python_packages.simple_physics import (
    GetContactBiasName,
    GetContactNodeModelName,
    contactcharge_edge,
)

from . import physics_t

# 과제 28쪽 유효 일함수 표 [eV]. 채점 고정값이므로 수정하지 않는다.
GATE_METALS = {
    "n+poly": 4.05,
    "Al": 4.10,
    "Ta": 4.25,
    "Ti": 4.33,
    "TaN": 4.45,
    "W": 4.60,
    "TiN": 4.65,
    "Mo": 4.70,
    "Ni": 5.10,
    "p+poly": 5.15,
    "Pt": 5.30,
}


def work_function(metal: str) -> float:
    """금속 이름 -> 유효 일함수 [eV]. 이름이 유일한 진실 원천이다."""
    try:
        return GATE_METALS[metal]
    except KeyError:
        raise ValueError(
            f"게이트 금속 '{metal}' 는 과제 28쪽 표에 없다. 허용: {sorted(GATE_METALS)}"
        ) from None


def gate_offset(metal: str, t_k: float, calibrate: bool = True) -> float:
    """게이트 접점 방정식에 더할 오프셋 [V] = phi_m - phi_i(T).

    phi_i 는 온도에 따라 내려가므로 (300 K 4.6137 eV -> 398 K 4.6007 eV)
    오프셋도 온도마다 다시 계산해야 한다. 300 K 값을 고정해 쓰면 고온 Vth 가
    너무 낮게, 고온 Ioff 가 너무 높게 나온다 (보수적이지만 틀린 값이다).
    """
    return work_function(metal) - physics_t.intrinsic_work_function(t_k, calibrate)


def CreateOxideContactWithWorkFunction(device, region, contact, offset_v):
    """simple_physics.CreateOxideContact 의 대체본. offset_v = phi_m - phi_i [V].

    offset_v = 0 이면 배포 코드와 완전히 같은 식이 되어 결과가 일치한다.
    """
    contact_bias_name = GetContactBiasName(contact)
    contact_model_name = GetContactNodeModelName(contact)
    eq = "Potential - {0} + ({1:+.12e})".format(contact_bias_name, float(offset_v))

    CreateContactNodeModel(device, contact, contact_model_name, eq)
    # 오프셋은 상수이므로 Potential 에 대한 미분은 배포 코드와 같이 1 이다.
    CreateContactNodeModelDerivative(device, contact, contact_model_name, eq, "Potential")

    if not InEdgeModelList(device, region, contactcharge_edge):
        CreateEdgeModel(device, region, contactcharge_edge, "Permittivity*ElectricField")
        CreateEdgeModelDerivatives(
            device, region, contactcharge_edge, "Permittivity*ElectricField", "Potential"
        )

    contact_equation(
        device=device,
        contact=contact,
        name="PotentialEquation",
        node_model=contact_model_name,
        edge_charge_model=contactcharge_edge,
    )
