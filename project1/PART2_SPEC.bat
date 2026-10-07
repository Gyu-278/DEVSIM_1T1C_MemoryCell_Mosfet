@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo ============================================================
echo  [P2-2] 1T1C 셀 스펙 측정 (개정: 0.5 ns 과도 읽기 + 리텐션)
echo ============================================================
echo  제출 구조(p+ 탭) 측정 + 바닥 접점 결과와 비교
echo  결과: results\bodytap\p2_G2_bodytap.json  ·  약 50분 (과도 해석 수백 회)
echo.
"%~dp0..\.conda\python.exe" -u -m part2.body_tap
if errorlevel 1 (
  echo.
  echo [FAILED] 위 오류 메시지를 확인하세요.
)
echo.
pause
