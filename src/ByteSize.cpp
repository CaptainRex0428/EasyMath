#include "ByteSize.h"

double EM::ByteSizeTo(int bytesize, SystemSizeUnit unit)
{
	// 定义转换常量
	const double KB = 1024.0;
	const double MB = 1024.0 * KB;
	const double GB = 1024.0 * MB;
	const double TB = 1024.0 * GB;
	const double PB = 1024.0 * TB;
	const double EB = 1024.0 * PB;
	const double ZB = 1024.0 * EB;
	const double YB = 1024.0 * ZB;


	// 根据目标单位进行转换
	switch (unit)
	{
	case SizeUnitByte:
		return bytesize; // 直接返回字节数
	case SizeUnitKB:
		return bytesize / KB; // 转换为KB
	case SizeUnitMB:
		return bytesize / MB; // 转换为MB
	case SizeUnitGB:
		return bytesize / GB; // 转换为GB
	case SizeUnitTB:
		return bytesize / TB; // 转换为TB
	case SizeUnitPB:
		return bytesize / PB; // 转换为PB
	case SizeUnitEB:
		return bytesize / EB; // 转换为EB
	case SizeUnitZB:
		return bytesize / ZB; // 转换为ZB
	case SizeUnitYB:
		return bytesize / YB; // 转换为YB
	default:
		return bytesize; // 默认返回字节数
	}

}

std::tuple<const char*, double> EM::ByteSizeConvert(uintmax_t bytesize)
{
	double size = (double)bytesize;
	const char* unit = "Byte";

	// to KB
	double sizebuf = size / 1024;
	if (((int)sizebuf) == 0)
	{
		return { unit,size };
	}
	size = sizebuf;
	unit = "KB";

	// to MB
	sizebuf = size / 1024;
	if (((int)sizebuf) == 0)
	{
		return { unit,size };
	}
	size = sizebuf;
	unit = "MB";

	// to GB
	sizebuf = size / 1024;
	if (((int)sizebuf) == 0)
	{
		return { unit,size };
	}
	size = sizebuf;
	unit = "GB";

	// to TB
	sizebuf = size / 1024;
	if (((int)sizebuf) == 0)
	{
		return { unit,size };
	}
	size = sizebuf;
	unit = "TB";

	// to PB
	sizebuf = size / 1024;
	if (((int)sizebuf) == 0)
	{
		return { unit,size };
	}
	size = sizebuf;
	unit = "PB";

	// to EB
	sizebuf = size / 1024;
	if (((int)sizebuf) == 0)
	{
		return { unit,size };
	}
	size = sizebuf;
	unit = "EB";

	// to ZB
	sizebuf = size / 1024;
	if (((int)sizebuf) == 0)
	{
		return { unit,size };
	}
	size = sizebuf;
	unit = "PB";

	return { unit,size };
}
