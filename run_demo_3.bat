@echo off
REM ==========================================================
REM  Genisys2 - one-click demo: build, test, live demo, visuals
REM  Usage: double-click, or run  run_demo_3.bat  in the repo
REM ==========================================================
cd /d "%~dp0"

echo.
echo ===== 1/5  COMPILING =====
g++ -std=c++17 -O2 tests/distance_tests_2.cpp -o test_distance.exe
if errorlevel 1 goto fail
g++ -std=c++17 -O2 tests/tree_validation_3.cpp -o test_tree.exe
if errorlevel 1 goto fail
g++ -std=c++17 -O2 src/genisys_app_3.cpp -o genisys.exe
if errorlevel 1 goto fail
echo Compiled OK.

echo.
echo ===== 2/5  UNIT + EDGE-CASE TESTS =====
test_distance.exe
if errorlevel 1 goto fail
test_tree.exe
if errorlevel 1 goto fail

echo.
echo ===== 3/5  LIVE DEMO (200 reads, Levenshtein) =====
echo Menu choices used: 6=stats, 3=range search, 4=k-nearest (+reverse complement), 5=VP-tree vs brute force
genisys.exe data\sample_reads_3.fasta levenshtein < demo_input_3.txt

echo.
echo ===== 4/5  GIT COMMIT HISTORY (per member) =====
git shortlog -sn
git log --format="%h | %an | %s" -n 12

echo.
echo ===== 5/5  OPENING VISUALIZERS =====
start "" visual\vp_tree_explained_1.html
start "" visual\dna_tree_visualizer_1.html

echo.
echo ALL DONE - everything passed.
pause
exit /b 0

:fail
echo.
echo *** SOMETHING FAILED - read the message above ***
pause
exit /b 1
