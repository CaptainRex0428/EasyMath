@echo off
pushd ..
Scripts\Premake\premake5.exe --file=Build.lua vs2022
popd
pause