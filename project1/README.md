# Project #1 — 로컬 실행 안내 (Part 1 + Part 2)

## 놓을 위치

`project1\` 을 `.conda`, `step2` 가 있는 폴더에 둔다.

```
...\DEVSIM_PROJECT\DEVSIM_PROJECT\
├─ .conda\                    (이미 있음)
├─ step2\                     (이미 있음, 하나도 안 고쳤다)
└─ project1\
   ├─ part1\  part2\           코드
   ├─ check_structure_file\    공식 검사기
   ├─ PART1_*.bat  PART2_*.bat 더블클릭 실행
   ├─ results\                 이미 잰 결과 (CSV·그림). 제출 구조 결과는 results\bodytap\
   ├─ report\                  보고서 그림(fig\), 그림 데이터(data\), 그림 스크립트(tools\)
   ├─ part1_2021144048.devsim  제출물, 검사기 통과
   └─ part2_2021144048.devsim  제출물, 검사기 통과
```

`project1` 이 반드시 `.conda` 와 **같은 폴더**에 있어야 한다. bat 이 `..\.conda\python.exe` 를,
코드가 `..\step2` 를 찾는다.

## 실행: bat 더블클릭


| bat 파일               | 하는 일                           | 결과                                       | 걸리는 시간 |
| -------------------- | ------------------------------ | ---------------------------------------- | ------ |
| `PART1_EXPORT.bat`   | 제출 파일 저장(p+ 탭 구조) + 공식 검사기     | `part1_2021144048.devsim`, 화면에 `OK`      | 10초    |
| `PART1_SPEC.bat`     | 최종 설계 8항목(탭 구조) + 바닥 접점 결과와 비교 | `results\bodytap\p1_final_bodytap.json`  | 약 6분   |
| `PART1_FIGURES.bat`  | 구조 그림 · Id-Vg 곡선 (탭 구조)        | `results\fig_*.html`, `curves_final.csv` | 5분     |
| `PART1_PHYSICS.bat`  | 추가한 물리 전후 비교 (보고서 ①)           | `results\physics_check_baseline.csv`     | 10분    |
| `PART1_COLLECT.bat`  | 탐색 기록표 모으기 (보고서 ②)             | `results\search_record.csv`              | 즉시     |
| `PART1_ABLATION.bat` | 파라미터별 ablation (보고서 ③)         | `results\ablation_final.csv`             | 1시간    |




### Part 2


| bat 파일              | 하는 일                                 | 결과                                               | 걸리는 시간 |
| ------------------- | ------------------------------------ | ------------------------------------------------ | ------ |
| `PART2_EXPORT.bat`  | 제출 파일 저장(p+ 탭 구조) + 공식 검사기           | `part2_2021144048.devsim`, 화면에 `OK`              | 10초    |
| `PART2_SPEC.bat`    | 과도 읽기 + 리텐션 전 항목(탭 구조) + 바닥 접점 비교    | `results\bodytap\p2_G2_bodytap.json`             | 약 50분  |
| `PART2_FIGURES.bat` | 셀 구조 그림 · ΔV_BL(t)/I(t) 과도 곡선 (탭 구조) | `results\fig_p2_*.html`, `p2_read_curves_G2.csv` | 3분     |


`PART1_PHYSICS` · `PART1_COLLECT` · `PART1_ABLATION` 과 `part2.spec_p2 --cell <이름>` 은 탐색 기록이라 바닥 접점 구조로 잰다.



## 설계를 바꿔 보고 싶으면

`part1\designs.py` 의 숫자만 고치면 됩니다. 그 다음 `PART1_SPEC.bat` 을 다시 돌리세요.

```python
CAND_A = Part1Device(
    gate_length_um=0.50,        # Lg
    oxide_thickness_nm=5.5,     # tox
    junction_depth_um=0.05,     # xj
    silicon_thickness_um=0.50,  # t_Si
    body_doping_cm3=7.5e16,     # NA
    sd_doping_cm3=1.0e20,       # ND
    gate_material="TaN",        # 28쪽 표의 이름만 (숫자 아님)
)
```

- **t_Si 는 드레인 공핍층 폭의 2배 이상**이어야 함. 얇으면 공핍층이 기판 바닥 경계까지 닿아서 시뮬레이션이 수렴하지 않음.
- **Lg 0.3 µm 는 과제가 허용하는 최솟값이지만 실제로는 punch-through** 가 납니다. 0.5 µm 이상을 활용하기.



## 최종 설계



### Part 1 (스펙 소자)

Lg 0.50 µm · tox 5.0 nm · xj 0.05 µm · t_Si 0.50 µm · S/D 각 0.50 µm ·
NA 7.5e16 · ND 1e20 · 게이트 TaN (4.45 eV)

### Part 2 (DRAM 1T1C 셀, 개정 스펙)

TRPASS: Lg 0.30 µm(고정) · tox 5.1 nm · xj 0.10 µm · t_Si 1.0 µm · S/D 각 0.25 µm ·
NA 9.5e16 · ND 9e20 · 게이트 Mo (4.70 eV) · body 접점 p+ 탭 위 (폭 0.15 µm, 소스와 간격 0.30 µm)
CAP: ZrO2 두 장 · 높이 1.40 µm · 두께 4.5 nm · 기둥 폭 40 nm · W = 100 nm

(Part 1 도 p+ 탭 위 body 접점: 폭 0.15 µm, 소스와 간격 0.25 µm. 8항목은 탐색 구조와 넷째 자리까지 같고 125 °C 누설만 20.77 → 20.69 pA/µm.)


| 항목                    | 기준        | 결과                                           |
| --------------------- | --------- | -------------------------------------------- |
| 읽기 '0' (0.5 ns)       | ≥ 60 mV   | 132.6                                        |
| 읽기 '1' (0.5 ns)       | ≥ 60 mV   | 69.61                                        |
| 64 ms 후 읽기 '0'        | ≥ 60 mV   | 125.0                                        |
| 64 ms 후 읽기 '1'        | ≥ 60 mV   | 66.54                                        |
| 리텐션 '1'               | ≥ 64 ms   | 185.6                                        |
| 리텐션 '0'               | ≥ 64 ms   | V_crit 미도달 (10 s 까지 0.205 V, V_crit 0.514 V) |
| C_STORE (plate dQ/dV) | ≤ 20 fF   | 19.24                                        |
| E_ox (2.5 V)          | ≤ 5 MV/cm | 4.902                                        |
| E_diel (1 V)          | ≤ 4 MV/cm | 2.222                                        |




## 파일별 역할


| 파일                                                   | 내용                                                                                |
| ---------------------------------------------------- | --------------------------------------------------------------------------------- |
| `part1\physics_t.py`                                 | 추가 물리 1. Eg(T), Nc/Nv(T), n_i(T), 전위 기준 phi_i(T)                                  |
| `part1\gate_physics.py`                              | 추가 물리 2. 게이트 일함수가 들어간 산화막 접점                                                      |
| `part1\simulator.py`                                 | step2 시뮬레이터 상속. 위 둘 + gate_metal 영역 + 적응 램프 + 128비트 정밀도                           |
| `part1\device.py`                                    | 설계 변수                                                                             |
| `part1\designs.py`                                   | 후보 설계                                                                             |
| `part1\spec.py`                                      | 8항목 측정 규칙                                                                         |
| `part1\search.py` `refine.py` `collect.py`           | 설계 탐색                                                                             |
| `part1\ablation.py` `physics_check.py`               | 보고서용 비교표                                                                          |
| `part1\figures.py`                                   | 구조 그림 · 특성 곡선                                                                     |
| `part1\export_devsim.py`                             | 제출 파일 생성(탭 구조) + 검사기 실행                                                           |
| `part1\body_tap.py`                                  | 제출 구조: p+ 탭 위 body 접점(TapMixin, TapPart1Simulator), 최종 8항목 측정                     |
| `part2\cell.py`                                      | 1T1C 셀. 유전막 hk_l/hk_r, 기둥, storage/plate 접점, dQ/dV 용량. 탐색한 후보 전부(CELLS)와 최종 FINAL |
| `part2\transient.py`                                 | **개정 스펙** Δt 스텝 과도 읽기(10 ps) + 적응 스텝 리텐션 엔진                                       |
| `part2\spec_p2.py`                                   | **개정 스펙** 전 항목 측정                                                                 |
| `part2\figures_p2.py`                                | **개정 스펙** ΔV_BL(t)·I(t) 과도 곡선, 구조 그림                                              |
| `part2\spec2.py` `search2.py` `figures2.py`          | 9/23 구버전 스펙용. 개정으로 쓰지 않음(기록용으로 남김)                                                |
| `part2\ranges.py`                                    | 과제 허용 범위 검사. 공식 검사기가 안 보는 도핑·치수를 막는다                                              |
| `part2\body_tap.py`                                  | 제출 구조: p+ 탭 위 body 접점(TapPart2Simulator), 최종 읽기·리텐션 측정                            |
| `part2\search2.py` `figures2.py` `export_devsim2.py` | 측정 · 그림 · 제출 파일                                                                   |


`step2` 는 한 줄도 고치지 않았습니다. Part 1 에서 더한 물리는 전부 `part1` 안에 있습니다.