#include "ByteSize.h"

double EasyMath::ByteSizeTo(int bytesize, SystemSizeUnit unit)
{

	const double KB = 1024.0;
	const double MB = 1024.0 * KB;
	const double GB = 1024.0 * MB;
	const double TB = 1024.0 * GB;
	const double PB = 1024.0 * TB;
	const double EB = 1024.0 * PB;
	const double ZB = 1024.0 * EB;
	const double YB = 1024.0 * ZB;



	switch (unit)
	{
	case SizeUnitByte:
		return bytesize; // ֱ
	case SizeUnitKB:
		return bytesize / KB; //
	case SizeUnitMB:
		return bytesize / MB; // 
	case SizeUnitGB:
		return bytesize / GB; // 
	case SizeUnitTB:
		return bytesize / TB; // 
	case SizeUnitPB:
		return bytesize / PB; // 
	case SizeUnitEB:
		return bytesize / EB; // 
	case SizeUnitZB:
		return bytesize / ZB; // 
	case SizeUnitYB:
		return bytesize / YB; // 
	default:
		return bytesize; // 
	}

}

std::tuple<const char*, double> EasyMath::ByteSizeConvert(uintmax_t bytesize)
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
