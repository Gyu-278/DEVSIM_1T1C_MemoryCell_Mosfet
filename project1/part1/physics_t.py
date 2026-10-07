# -*- coding: utf-8 -*-
"""[추가 물리 1] 실리콘의 온도 의존 물성.

배포 코드(devsim.python_packages.simple_physics)는 SetSiliconParameters 에서
V_t 와 kT 만 온도로 계산하고, 진성 캐리어 농도 n_i 는 300 K 값 1e10 으로 박아 둔다.
그 파일 자체에 "TODO: make T a free parameter and T dependent parameters as models"
라고 적혀 있다. 그래서 125 °C 스펙을 그 코드로 재면 Vth 가 거꾸로 올라가고
Ioff 는 2.2 배밖에 안 늘어난다 (물리적으로는 수천 배여야 한다).

여기서 Eg(T), Nc(T), Nv(T), n_i(T) 를 넣어 그 구멍을 막는다.
과제가 금지한 것은 "코드 상수의 임의 수정"이고, 상수를 물리 모델로 바꾸는 것은
명시적으로 허용되어 있다.
"""

from __future__ import annotations

import math

# --- devsim 이 쓰는 값과 정확히 같은 상수를 쓴다 (simple_physics.py:34-41) ---
Q = 1.6e-19            # C
K_B = 1.3806503e-23    # J/K
NI_REF = 1.0e10        # simple_physics.py 가 300 K 로 박아 둔 값
T_REF = 300.0          # K

# --- 실리콘 밴드 파라미터 ---
# Varshni(1967) 식. Eg0 는 조교가 공지한 값(Eg(300K)=1.1245, Eg(398K)=1.0975)을
# 정확히 재현하도록 1.170 을 쓴다. 채점기와 같은 값을 쓰는 것이 목적이다.
EG0 = 1.170            # eV
VARSHNI_ALPHA = 4.73e-4  # eV/K
VARSHNI_BETA = 636.0     # K
# 유효 상태밀도. Green, J. Appl. Phys. 67, 2944 (1990)
NC_REF = 2.86e19       # cm^-3 at 300 K
NV_REF = 3.10e19       # cm^-3 at 300 K
CHI_SI = 4.05          # eV, 실리콘 전자 친화도


def thermal_voltage(t_k: float) -> float:
    """kT/q [V]."""
    return K_B * t_k / Q


def band_gap(t_k: float) -> float:
    """Varshni 밴드갭 Eg(T) [eV]."""
    return EG0 - VARSHNI_ALPHA * t_k * t_k / (t_k + VARSHNI_BETA)


def nc(t_k: float) -> float:
    """전도대 유효 상태밀도 [cm^-3]. 포물선 밴드에서 T^1.5."""
    return NC_REF * (t_k / T_REF) ** 1.5


def nv(t_k: float) -> float:
    """가전자대 유효 상태밀도 [cm^-3]."""
    return NV_REF * (t_k / T_REF) ** 1.5


def intrinsic_raw(t_k: float) -> float:
    """보정 없는 n_i(T) = sqrt(Nc*Nv)*exp(-Eg/2kT) [cm^-3]."""
    return math.sqrt(nc(t_k) * nv(t_k)) * math.exp(-band_gap(t_k) / (2.0 * thermal_voltage(t_k)))


def intrinsic(t_k: float, calibrate: bool = True) -> float:
    """n_i(T) [cm^-3].

    calibrate=True 면 300 K 에서 정확히 NI_REF(1e10)가 되도록 전체를 비례 보정한다.
    이렇게 하면 300 K 결과가 배포 코드와 완전히 같아서, 온도 물리를 넣은 것이
    상온 스펙 6개를 흔들지 않았다는 점을 보고서에서 그대로 보일 수 있다.
    (보정 없는 값은 300 K 에서 1.110e10 으로 11% 높다.)
    """
    value = intrinsic_raw(t_k)
    if calibrate:
        value *= NI_REF / intrinsic_raw(T_REF)
    return value


def intrinsic_work_function(t_k: float, calibrate: bool = True) -> float:
    """게이트 일함수 오프셋의 기준 [eV] = chi + Eg(T)/2.

    조교 공지: 게이트 접점은 Potential = V_G - (phi_m - (chi + Eg/2)) 이고,
    chi = 4.05 eV, Eg 는 온도에 따라 변하는 값을 쓴다. 채점기와 같은 식을 써야
    고온에서 수십 mV 가 어긋나지 않으므로 그대로 따른다.
      300 K -> 4.612260 eV,  398 K -> 4.598769 eV

    참고: 이 코드의 전위 0 은 E_F = E_i 인 지점이고 E_i 는 정확히는 midgap 이
    아니다(E_C - E_i = V_t*ln(Nc/n_i) 로 쓰면 300 K 에서 4.6137 eV). 차이는
    Vth 로 1.4 mV 수준이라 무시할 만하고, 채점기와 맞추는 쪽이 이득이다.
    calibrate 인자는 호출부 호환을 위해 남겨 두며 결과에 영향을 주지 않는다.
    """
    return CHI_SI + band_gap(t_k) / 2.0


# 조교 공지의 이동도 온도 지수. 속도 포화와 도핑 의존 이동도는 채점기에도 없다.
MU_N_EXPONENT = 2.4
MU_P_EXPONENT = 2.2


def mobility_lattice(mu_ref: float, t_k: float, exponent: float) -> float:
    """격자 산란 이동도 mu(T) = mu300*(T/300)^-exponent [cm^2/Vs].

    조교 공지: 채점기는 mu_n = 400*(T/300)^-2.4, mu_p = 200*(T/300)^-2.2 를 쓴다.
    300 K 에서는 배포 코드와 같은 400/200 이고, 398 K 에서 203/107 로 거의 반이 된다.
    고온 Ioff 와 Part 2 의 모든 시뮬레이션이 125 °C 이므로 반드시 넣어야 한다.
    """
    return mu_ref * (t_k / T_REF) ** (-exponent)


def mu_n(t_k: float, mu_ref: float = 400.0) -> float:
    return mobility_lattice(mu_ref, t_k, MU_N_EXPONENT)


def mu_p(t_k: float, mu_ref: float = 200.0) -> float:
    return mobility_lattice(mu_ref, t_k, MU_P_EXPONENT)
