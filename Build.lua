include "Directory.lua"
include "Dependencies.lua"
include "Config.lua"

workspace "EasyMath"
	architecture "x64"
	startproject "Sandbox"
	configurations
	{
		"Debug",
		"Release",
        "Dist"
	}

    filter "system:windows"
    buildoptions { "/EHsc", "/Zc:preprocessor", "/Zc:__cplusplus","/utf-8" }
	
	add_rider_refreshconfig()

group ""
	include "EasyMath.lua"

group "Sandbox"
	include "Sandbox.lua"