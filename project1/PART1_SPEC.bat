@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo ============================================================
echo  [1] 최종 설계 8항목 스펙 측정
echo ============================================================
echo  제출 구조(p+ 탭) 8항목 + 바닥 접점 결과와 비교 + 검사기
echo  결과: results\bodytap\p1_final_bodytap.json  ·  약 6분
echo.
"%~dp0..\.conda\python.exe" -u -m part1.body_tap
if errorlevel 1 (
  echo.
  echo [FAILED] 위 오류 메시지를 확인하세요.
)
echo.
pause
