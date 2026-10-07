# -*- coding: utf-8 -*-
"""[Project#1 Part 1] HW#1(step2) 시뮬레이터에 빠진 물리를 추가한 패키지.

step2 는 완료된 과제라 손대지 않는다. 여기서는 step2 의 MosfetSimulator 를
상속해서 Part 1 에 필요한 물리만 얹는다.
"""

import os
import sys
from pathlib import Path

# 이 문제 크기에서는 MKL/OpenMP 멀티스레딩이 순손실이다. 측정값: 기본 5.63 s에
# CPU 142%, 스레드 1개로 제한하면 4.92 s에 CPU 22%. 스레드 조율 비용이 이득보다 크다.
# 변형을 여러 프로세스로 병렬 측정할 때는 이 차이가 6배 이상의 처리량으로 돌아온다.
# numpy/devsim 이 import 시점에 읽으므로 반드시 그 전에 설정해야 한다.
# 굳이 스레드를 쓰고 싶으면 DEVSIM_KEEP_THREADS=1 을 걸면 된다.
if not os.environ.get("DEVSIM_KEEP_THREADS"):
    for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(_v, "1")

# step2 를 import 경로에 추가한다 (project1/part1 -> DEVSIM_PROJECT/step2)
_STEP2 = Path(__file__).resolve().parent.parent.parent / "step2"
if str(_STEP2) not in sys.path:
    sys.path.insert(0, str(_STEP2))
