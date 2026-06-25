@echo off
REM Batch extract all RTR4 chapters to _raw/
REM Run from E:\Private\EasyMath

setlocal enabledelayedexpansion

set "RAW=E:\Private\EasyMath\.claude\memory\reference\RTR4\_raw"
set "SRC=E:\Private\EasyMath\reference\Real-Time Rendering, Fourth Edition (Tomas Akenine-Moller, Eric Haines, Naty Hoffman etc.)\EPUB\text"
set "SCRIPT=E:\Private\EasyMath\.claude\memory\reference\_scripts\extract_epub_chapter.py"

for /L %%i in (1,1,26) do (
    set "NUM=00%%i"
    set "NUM=!NUM:~-3!"
    if exist "%SRC%\ch!NUM!.xhtml" (
        set "EXT=xhtml"
    ) else (
        set "EXT=html"
    )
    set "INFILE=%SRC%\ch!NUM!.!EXT!"
    set "OUTFILE=%RAW%\ch!NUM!.txt"
    echo Extracting ch!NUM!.!EXT! ...
    py "%SCRIPT%" "!INFILE!" "!OUTFILE!"
)

echo.
echo All chapters extracted to %RAW%
dir "%RAW%"
