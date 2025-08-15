LocationDir = "%{wks.location}/settings"

DependenciesRoot = "%{wks.location}/Extension/MultiExtend/Dependencies"
ExtensionRoot = "%{wks.location}/Extension"

OutputDir = "%{cfg.system}-%{cfg.architecture}/%{cfg.buildcfg}"

EXEDir = "%{wks.location}/bin/"..OutputDir
LibDir = "%{wks.location}/bin/"..OutputDir
DynamicDir = "%{wks.location}/bin/"..OutputDir
TargetDir = "%{wks.location}/bin/" .. OutputDir .. "/%{prj.name}"
ObjectDir = "%{wks.location}/bin/obj/" .. OutputDir .. "/%{prj.name}"