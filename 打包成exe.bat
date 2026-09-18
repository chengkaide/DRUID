@echo off
rem =====================================================================
rem  Build a standalone .exe
rem
rem  All comments are ASCII-only on purpose: cmd.exe decodes .bat files
rem  with the active code page, so non-ASCII bytes can turn to garbage.
rem
rem  Result:  dist\DRUID.exe   (~150-250 MB, single file)
rem           Double-click -> local web service starts -> browser opens.
rem =====================================================================
setlocal

set "PYEXE="
rem the virtualenv may be named "druid" (new) or "upb" (historical); accept both
for /d %%p in ("%USERPROFILE%\.workbuddy\binaries\python\envs\druid") do (
    set "PYEXE=%%~p\Scripts\python.exe"
)
if not exist "%PYEXE%" (
    for /d %%p in ("%USERPROFILE%\.workbuddy\binaries\python\envs\upb") do (
        set "PYEXE=%%~p\Scripts\python.exe"
    )
)
if not exist "%PYEXE%" (
    for /d %%u in ("%SYSTEMDRIVE%\Users\*") do (
        for %%e in (druid upb) do (
            for /d %%p in ("%%~u\.workbuddy\binaries\python\envs\%%e") do (
                if exist "%%~p\Scripts\python.exe" set "PYEXE=%%~p\Scripts\python.exe"
            )
        )
    )
)
if not exist "%PYEXE%" (
    echo [ERROR] Python not found.
    pause & exit /b 1
)

cd /d "%~dp0"
echo Building with: %PYEXE%
echo This takes several minutes. Please wait...
echo.

rem --collect-data druid does NOT work here: druid is not a pip-installed package,
rem so PyInstaller cannot locate its package data. The static assets have to be
rem added explicitly (Windows separator between source and dest is ';').
"%PYEXE%" -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onefile ^
  --console ^
  --name DRUID ^
  --paths . ^
  --add-data "druid/webui/static;druid/webui/static" ^
  --hidden-import=pandas ^
  --hidden-import=numpy ^
  --hidden-import=matplotlib ^
  --hidden-import=matplotlib.backends.backend_agg ^
  --hidden-import=openpyxl ^
  --hidden-import=xlrd ^
  --exclude-module=tkinter ^
  --exclude-module=IPython ^
  --exclude-module=PyQt5 ^
  --exclude-module=PySide2 ^
  --distpath dist ^
  --workpath build_tmp ^
  launcher.py

echo.
if exist "dist\DRUID.exe" (
    echo   [OK] dist\DRUID.exe
    for %%F in ("dist\DRUID.exe") do echo        size: %%~zF bytes
) else (
    echo   [FAILED] dist\DRUID.exe not found.
)
echo.
pause
