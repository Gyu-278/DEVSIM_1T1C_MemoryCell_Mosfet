#!/bin/bash
# 압축본 안에 project1/ 폴더째 들어간다. 로컬의 DEVSIM_PROJECT\DEVSIM_PROJECT\ 에서
# "여기에 풀기" 하면 project1\ 이 생기고 bat 들이 바로 동작한다 (..\.conda, ..\step2 를 쓴다).
cd "$(dirname "$0")/.."
find project1/part1 project1/part2 -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null
Z=project1_final.zip; rm -f $Z
zip -qr $Z project1/part1 project1/part2 project1/check_structure_file \
   project1/PART1_*.bat project1/PART2_*.bat project1/README.md project1/REPORT_CODE.md project1/make_zip.sh \
   project1/part1_2021144048.devsim project1/part2_2021144048.devsim \
   project1/results/*.csv project1/results/*.html project1/results/*.json project1/results/*.txt project1/results/refine \
   project1/results/bodytap \
   project1/report/Project1_보고서.docx project1/report/Project1_보고서.html project1/report/Project1_보고서.md project1/report/fig project1/report/tools \
   -x '*__pycache__*' '*.pyc'
echo "$(pwd)/$Z  $(du -h $Z | cut -f1)"
