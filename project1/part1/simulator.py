# -*- coding: utf-8 -*-
"""Part 1 시뮬레이터. step2 MosfetSimulator 를 상속해 네 가지를 더한다.

  1) 게이트 일함수      gate_physics.CreateOxideContactWithWorkFunction
  2) 온도 의존 n_i      physics_t.intrinsic  (n_i · SRH 의 n1 · p1 을 함께 갱신)
  3) gate_metal 영역    제출용 .devsim 파일이 요구하는 라벨 영역
  4) 적응 바이어스 램프  한 스텝이 안 풀리면 직전 점으로 되돌리고 스텝을 반으로 줄인다

RAMP_STEP_V · MU_N · MU_P 등 배포 코드의 상수는 바꾸지 않는다. 4)는 상수를 고치는
것이 아니라 실패한 스텝에서만 더 잘게 나누어 가는 경로 변경이다.
"""

from __future__ import annotations

import math

import devsim

# 일부 PC 에서 MKL 경로가 Newton 을 2~4 회 만에 "RelError 0" 으로 끝내고 Id 가
# 모든 Vg 에서 같은 값으로 나오는 증상이 보고되었다(조교 안내). umfshim 을
# import 하면 선형 솔버가 MKL Pardiso 에서 UMFPACK(custom)으로 바뀌어 그 경로를
# 피한다. turtle 에서 두 솔버를 대조한 결과 Vth·SS·Ion·Ioff 가 출력한 모든 자릿수
# 까지 같았다(Vth 0.451848, Ioff 0.0576599 pA/um). 속도도 비슷하다(77 s 대 73 s).
try:
    from devsim.umfpack import umfshim  # noqa: F401
except Exception:
    pass

from devsim.python_packages.model_create import CreateSolution
from devsim.python_packages.simple_physics import (
    CreateOxidePotentialOnly,
    CreateSiliconOxideInterface,
    CreateSiliconPotentialOnly,
    CreateSiliconPotentialOnlyContact,
    GetContactBiasName,
    SetOxideParameters,
    SetSiliconParameters,
)
from mosfet_tool.simulator import NM, UM, MosfetSimulator

from . import gate_physics, physics_t
from .device import Part1Device


class Part1Simulator(MosfetSimulator):
    """Part 1 용 2-D nMOS. 호출 순서는 step2 와 같다."""

    MAX_BISECTIONS = 6          # 적응 램프 재분할 한도
    DY_METAL = 25.0 * NM        # gate_metal 라벨 영역의 메시 간격

    @staticmethod
    def enable_extended_precision(on: bool = True) -> None:
        """DEVSIM 128비트 확장 정밀도. 과제 36쪽이 수렴 실패 시 쓰라고 안내한 방법이다.

        Vg = 0 에서 드레인을 2 V 까지 올릴 때 기본 정밀도로는 Newton 이
        RelError 3.9e-5 근처에서 한계 순환에 빠져 relative_error 1e-6 에 도달하지
        못한다. 확장 정밀도는 물리도 솔버 허용오차도 건드리지 않고 그 문제만 없앤다.
        """
        for name in ("extended_solver", "extended_model", "extended_equation"):
            devsim.set_parameter(name=name, value=bool(on))

    def __init__(self, device: Part1Device, name: str = "part1",
                 structure_file: str | None = None, calibrate_ni: bool = True,
                 extended_precision: bool = True,
                 use_ni_temperature: bool = True, use_work_function: bool = True,
                 use_mobility_t: bool = True):
        super().__init__(device, name)
        self.extended_precision = extended_precision
        # 보고서 ① 용 스위치. False 로 두면 배포 코드와 똑같이 동작한다.
        self.use_ni_temperature = use_ni_temperature
        self.use_work_function = use_work_function
        self.use_mobility_t = use_mobility_t
        self.dev: Part1Device = device
        self.y_metal_top = self.y_oxide_top - device.gate_metal_thickness_nm * NM
        self.structure_file = structure_file
        self.calibrate_ni = calibrate_ni
        self.ramp_retries = 0

    # ------------------------------------------------------------------
    # 1) 구조
    # ------------------------------------------------------------------
    def build(self) -> None:
        """메시·도핑까지 만든 직후에 구조 파일을 저장하고, 그 다음 물리를 얹는다.

        저장 시점이 중요하다. 물리를 등록한 뒤에 저장하면 방정식과 Potential 노드
        모델이 파일에 같이 들어가서 제출 검사기가 REJECT 한다.
        """
        self._clear_session()
        if self.extended_precision:
            self.enable_extended_precision(True)
        self._build_mesh()
        self._build_doping()
        if self.structure_file:
            self.save_structure(self.structure_file)
        self._build_physics()

    def save_structure(self, path: str) -> str:
        """구조만 담긴 .devsim 파일을 쓴다. 물리 등록 전에만 호출해야 한다."""
        devsim.write_devices(file=str(path), device=self.name, type="devsim")
        return str(path)

    def _build_mesh(self) -> None:
        """step2 메시에 gate_metal 라벨 영역을 위로 덧붙인다.

        메시 위쪽 여유를 먼저 넓히지 않으면 devsim 이 그 영역을 조용히 버린다.
        """
        pad = 1.0e-8
        x_min, x_max = -pad, self.x_right + pad
        y_min, y_max = self.y_metal_top - pad, self.y_bottom + pad

        devsim.create_2d_mesh(mesh=self.mesh)
        for pos in (x_min, 0.0, self.x_gate_left, self.x_gate_right, self.x_right, x_max):
            devsim.add_2d_mesh_line(mesh=self.mesh, dir="x", pos=pos, ps=self.DX_CHANNEL)
        for pos, spacing in (
            (y_min, self.DY_METAL),
            (self.y_metal_top, self.DY_METAL),
            (self.y_oxide_top, self.DY_OXIDE),
            (0.0, min(self.DY_OXIDE, self.DY_JUNCTION)),
            (self.y_junction, self.DY_JUNCTION),
            (self.y_bottom, self.DY_BULK),
            (y_max, self.DY_BULK),
        ):
            devsim.add_2d_mesh_line(mesh=self.mesh, dir="y", pos=pos, ps=spacing)

        # 순서 주의: 경계 없는 air 를 먼저 두고, 뒤의 영역이 그 위를 덮어쓴다.
        devsim.add_2d_region(mesh=self.mesh, material="Air", region="air")
        devsim.add_2d_region(mesh=self.mesh, material="Silicon", region="bulk",
                             xl=0.0, xh=self.x_right, yl=self.y_bottom, yh=0.0)
        devsim.add_2d_region(mesh=self.mesh, material="Oxide", region="oxide",
                             xl=self.x_gate_left, xh=self.x_gate_right,
                             yl=0.0, yh=self.y_oxide_top)
        # 게이트 재료 이름을 조교에게 전달하는 라벨 영역. 방정식은 절대 얹지 않는다.
        devsim.add_2d_region(mesh=self.mesh, material=self.dev.gate_material,
                             region="gate_metal",
                             xl=self.x_gate_left, xh=self.x_gate_right,
                             yl=self.y_oxide_top, yh=self.y_metal_top)

        # 게이트 접점은 규칙대로 oxide 영역 위에 둔다 (검사기: gate -> oxide).
        devsim.add_2d_contact(mesh=self.mesh, name="gate", region="oxide", material="metal",
                              xl=self.x_gate_left, xh=self.x_gate_right,
                              yl=self.y_oxide_top, yh=self.y_oxide_top)
        devsim.add_2d_contact(mesh=self.mesh, name="source", region="bulk", material="metal",
                              xl=x_min, xh=self.x_gate_left, yl=0.0, yh=0.0)
        devsim.add_2d_contact(mesh=self.mesh, name="drain", region="bulk", material="metal",
                              xl=self.x_gate_right, xh=x_max, yl=0.0, yh=0.0)
        devsim.add_2d_contact(mesh=self.mesh, name="body", region="bulk", material="metal",
                              xl=x_min, xh=x_max, yl=self.y_bottom, yh=self.y_bottom)
        # 계면은 방정식을 푸는 영역 사이에만 둔다. gate_metal 쪽에는 만들지 않는다.
        devsim.add_2d_interface(mesh=self.mesh, name="bulk_oxide",
                                region0="bulk", region1="oxide")
        devsim.finalize_mesh(mesh=self.mesh)
        devsim.create_device(mesh=self.mesh, device=self.name)

    # ------------------------------------------------------------------
    # 2) 물리
    # ------------------------------------------------------------------
    def _build_physics(self) -> None:
        t_k = self.dev.temperature_k
        # gate_metal 과 air 에는 Potential 을 만들지 않는다 (방정식 없는 라벨 영역).
        for region in ("bulk", "oxide"):
            CreateSolution(self.name, region, "Potential")

        SetSiliconParameters(self.name, "bulk", t_k)
        # [추가 물리 1] SetSiliconParameters 는 n_i 를 300 K 상수로 박아 둔다.
        # 모든 모델이 n_i 를 영역 파라미터로 읽으므로 파라미터만 덮어써도 충분하다.
        # SRH 의 n1 · p1 도 같은 값을 쓰므로 함께 갱신한다 (midgap trap 가정).
        if self.use_ni_temperature:
            self.n_i = physics_t.intrinsic(t_k, self.calibrate_ni)
            for pname in ("n_i", "n1", "p1"):
                devsim.set_parameter(device=self.name, region="bulk",
                                     name=pname, value=self.n_i)
        else:
            self.n_i = physics_t.NI_REF     # 배포 코드 그대로 (300 K 상수로 동결)
        # [추가 물리 3] 이동도의 온도 의존. 조교 공지의 채점기 식과 같다.
        #   mu_n = 400*(T/300)^-2.4,  mu_p = 200*(T/300)^-2.2
        # 300 K 에서는 배포 코드와 같은 400/200 이라 상온 스펙은 흔들리지 않고,
        # 398 K 에서 203/107 로 떨어져 고온 Ioff 와 Part 2 전체에 영향을 준다.
        # 상수를 임의로 바꾼 것이 아니라 상수를 물리 모델로 대체한 것이고,
        # 과제가 그 대체를 명시적으로 허용한다.
        self.mu_n = physics_t.mu_n(t_k, self.MU_N) if self.use_mobility_t else self.MU_N
        self.mu_p = physics_t.mu_p(t_k, self.MU_P) if self.use_mobility_t else self.MU_P
        devsim.set_parameter(device=self.name, region="bulk", name="mu_n", value=self.mu_n)
        devsim.set_parameter(device=self.name, region="bulk", name="mu_p", value=self.mu_p)

        CreateSiliconPotentialOnly(self.name, "bulk")
        SetOxideParameters(self.name, "oxide", t_k)
        CreateOxidePotentialOnly(self.name, "oxide", "log_damp")

        # [추가 물리 2] 일함수 오프셋이 들어간 게이트 접점. 온도마다 다시 계산한다.
        self.gate_offset_v = gate_physics.gate_offset(
            self.dev.gate_material, t_k, self.calibrate_ni) if self.use_work_function else 0.0
        gate_physics.CreateOxideContactWithWorkFunction(
            self.name, "oxide", "gate", self.gate_offset_v)
        devsim.set_parameter(device=self.name, name=GetContactBiasName("gate"), value=0.0)

        for contact in ("source", "drain", "body"):
            CreateSiliconPotentialOnlyContact(self.name, "bulk", contact)
            devsim.set_parameter(device=self.name, name=GetContactBiasName(contact), value=0.0)
        CreateSiliconOxideInterface(self.name, "bulk_oxide")

    # ------------------------------------------------------------------
    # 3) 적응 바이어스 램프
    # ------------------------------------------------------------------
    def set_bias(self, contact: str, volts: float) -> None:
        """RAMP_STEP_V 간격으로 올리되, 안 풀리는 스텝은 반으로 줄여 다시 시도한다."""
        self._ramp(contact, float(volts), self.RAMP_STEP_V, 0)

    def _ramp(self, contact: str, target: float, step: float, depth: int) -> None:
        start = self.bias[contact]
        if abs(target - start) < 1.0e-12:
            return
        n = max(1, math.ceil(abs(target - start) / step - 1.0e-9))
        for i in range(1, n + 1):
            value = start + (target - start) * i / n
            devsim.set_parameter(device=self.name,
                                 name=GetContactBiasName(contact), value=value)
            try:
                self._solve()
            except Exception:
                if depth >= self.MAX_BISECTIONS:
                    raise
                # 직전 수렴점으로 되돌린 뒤 더 잘게 나누어 다시 올라간다.
                back = start + (target - start) * (i - 1) / n
                devsim.set_parameter(device=self.name,
                                     name=GetContactBiasName(contact), value=back)
                self._solve()
                self.bias[contact] = back
                self.ramp_retries += 1
                self._ramp(contact, target, step / 2.0, depth + 1)
                return
            self.bias[contact] = value
        self.bias[contact] = target
