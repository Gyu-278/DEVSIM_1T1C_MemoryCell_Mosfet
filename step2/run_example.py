# -*- coding: utf-8 -*-
"""시뮬레이터 클래스를 직접 불러와 쓰는 최소 예시.

mosfet.py(CLI)를 거치지 않고 MosfetSimulator를 바로 다룬다.
클래스 호출 순서를 눈으로 확인하고 싶을 때 실행해 볼 것.
"""

from mosfet_tool.config import Device
from mosfet_tool.simulator import MosfetSimulator

device = Device()  # 기본값 사용. Device(gate_length_um=0.5)처럼 바꿀 수 있다.

sim = MosfetSimulator(device, name="example")
sim.build()
sim.solve_equilibrium()
sim.enable_transport()

curve = sim.sweep_idvg(vd=0.05, start=0.0, stop=2.0, step=0.1)  # pandas DataFrame

print("\n== Id-Vg (Vd = 0.05 V) ==")
print(curve.to_string(index=False))
