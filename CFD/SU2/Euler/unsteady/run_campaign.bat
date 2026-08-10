@echo off
rem ============================================================================
rem  run_campaign.bat -- the whole SU2 prescribed-pitch campaign, one command.
rem  Run from CFD\SU2\Euler\unsteady with the mflco env active:
rem      run_campaign.bat
rem  Steps: mesh check -> steady trim (AOA 10) -> pitch 2 deg -> pitch 14 deg
rem         -> post-process (F8 + cycle work). Everything logged to *.log so
rem  the run can be left overnight; rerun the script to redo only what failed
rem  (each step is skipped if its product already exists).
rem ============================================================================
rem run from anywhere: work inside this script's own folder
pushd "%~dp0"

rem --- 0. mesh must sit beside the cfgs --------------------------------------
if not exist naca0020_medium.su2 (
    if exist ..\naca0020_medium.su2 (
        copy ..\naca0020_medium.su2 . >nul
        echo [mesh] copied naca0020_medium.su2 from parent folder
    ) else (
        echo [ERROR] naca0020_medium.su2 not found here or in parent -- aborting.
        popd
        exit /b 1
    )
)

rem --- 1. steady trim solution (minutes) --------------------------------------
if exist restart_steady10.dat (
    echo [1/4] steady_a10: restart_steady10.dat exists -- skipping
) else (
    echo [1/4] steady_a10 running LIVE below -- watch CL settle near 1.1-1.2
    echo        (an interrupted earlier attempt? delete restart_steady10.dat first)
    SU2_CFD steady_a10.cfg
    if not exist restart_steady10.dat (
        echo [ERROR] steady run produced no restart file
        popd
        exit /b 1
    )
)
copy /y restart_steady10.dat solution_flow.dat >nul

rem --- 2. small amplitude: attached benchmark (1-4 h) -------------------------
if exist history_small.csv (
    echo [2/4] pitch_small: history_small.csv exists -- skipping
) else (
    echo [2/4] pitch_small_k0286 running -- SILENT BY DESIGN, ~1-4 h.
    echo        live view from a second cmd window:
    echo          powershell Get-Content pitch_small.log -Wait -Tail 5
    SU2_CFD pitch_small_k0286.cfg > pitch_small.log 2>&1
    if not exist history_small.csv (
        echo [WARN] no history_small.csv -- check pitch_small.log tail:
        powershell -NoProfile -command "Get-Content pitch_small.log -Tail 5"
    )
)

rem --- 3. large amplitude: ceiling demonstration (1-4 h) ----------------------
if exist history_large.csv (
    echo [3/4] pitch_large: history_large.csv exists -- skipping
) else (
    echo [3/4] pitch_large_k0286 running -- SILENT BY DESIGN, ~1-4 h.
    echo        live view: powershell Get-Content pitch_large.log -Wait -Tail 5
    SU2_CFD pitch_large_k0286.cfg > pitch_large.log 2>&1
    if not exist history_large.csv (
        echo [WARN] no history_large.csv -- check pitch_large.log tail:
        powershell -NoProfile -command "Get-Content pitch_large.log -Tail 5"
    )
)

rem --- 4. post-process: F8 + cycle work ---------------------------------------
echo [4/4] post-processing...
set PYTHONPATH=%~dp0..\..\..\..\src
python su2_unsteady_energy.py
if errorlevel 1 echo [WARN] post-processing failed -- run it from an mflco cmd window
echo.
echo Campaign done. Products: F8_euler_loops.png, su2_unsteady_energy.json
echo Logs: steady_a10.log, pitch_small.log, pitch_large.log
popd