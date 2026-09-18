@echo off
rem =====================================================================
rem  DRUID -- LA-ICP-MS zircon U-Pb reduction tool  (web launcher)
rem
rem  Usage:  double-click this file.  A browser page will open.
rem          Close this black window to quit.
rem
rem  NOTE: All Chinese text is deliberately avoided in this .bat file.
rem        Windows cmd.exe parses .bat files using the current code page,
rem        so non-ASCII characters may turn into garbage depending on how
rem        this file was saved. Every Chinese-looking path below is reached
rem        through variables that are expanded at RUN time (e.g.
rem        %USERPROFILE%), never written literally -- that keeps this file
rem        100% ASCII and therefore safe under any code page.
rem =====================================================================
setlocal enabledelayedexpansion
title DRUID - zircon U-Pb reduction

rem ---- 1. locate the managed Python interpreter -----------------------
rem Std-path variant first (fast path). %USERPROFILE% is expanded at run
rem time, so a Chinese user folder name causes no encoding trouble here.
rem
rem NOTE: the virtualenv directory may be named either "druid" (if it was
rem       created after the tool was renamed) or "upb" (its historical
rem       name). Both are accepted so that renaming the tool does not
rem       invalidate an existing environment.
set "PYEXE="
for /d %%p in ("%USERPROFILE%\.workbuddy\binaries\python\envs\druid") do (
    set "PYEXE=%%~p\Scripts\python.exe"
)
if not exist "%PYEXE%" (
    for /d %%p in ("%USERPROFILE%\.workbuddy\binaries\python\envs\upb") do (
        set "PYEXE=%%~p\Scripts\python.exe"
    )
)

rem Fallback: sweep every user profile. Useful when the file is copied to
rem a machine whose account name differs from the one it was made on.
if not exist "%PYEXE%" (
    for /d %%u in ("%SYSTEMDRIVE%\Users\*") do (
        for %%e in (druid upb) do (
            for /d %%p in ("%%~u\.workbuddy\binaries\python\envs\%%e") do (
                if exist "%%~p\Scripts\python.exe" set "PYEXE=%%~p\Scripts\python.exe"
            )
        )
    )
)

rem Last resort: whatever "python" resolves to on PATH.
if not exist "%PYEXE%" (
    for /f "delims=" %%c in ('where python 2^>nul') do (
        if not defined PYEXE set "PYEXE=%%c"
    )
)

if not exist "%PYEXE%" if not defined PYEXE (
    echo.
    echo   [ERROR] Python interpreter not found.
    echo.
    echo   Looked for the managed environment under:
    echo     %%USERPROFILE%%\.workbuddy\binaries\python\envs\druid
    echo     %%USERPROFILE%%\.workbuddy\binaries\python\envs\upb
    echo.
    pause
    exit /b 1
)

echo   Python: %PYEXE%

rem ---- 2. move into the package root so "python -m druid..." resolves ----
rem %~dp0 is this file's directory (contains the  druid/  package).
cd /d "%~dp0"

rem ---- 3. start the local web service ---------------------------------
echo.
echo   Starting local server... your browser will open in a moment.
echo   Keep this window open while you work. Close it to stop the service.
echo.
"%PYEXE%" -m druid.cli.serve_ui --port 0 --no-browser 2>&1

echo.
echo   Service stopped.
pause
