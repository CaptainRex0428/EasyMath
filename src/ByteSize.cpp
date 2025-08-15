#include "ByteSize.h"

double EM::ByteSizeTo(int bytesize, SystemSizeUnit unit)
{
	// ����ת������
	const double KB = 1024.0;
	const double MB = 1024.0 * KB;
	const double GB = 1024.0 * MB;
	const double TB = 1024.0 * GB;
	const double PB = 1024.0 * TB;
	const double EB = 1024.0 * PB;
	const double ZB = 1024.0 * EB;
	const double YB = 1024.0 * ZB;


	// ����Ŀ�굥λ����ת��
	switch (unit)
	{
	case SizeUnitByte:
		return bytesize; // ֱ�ӷ����ֽ���
	case SizeUnitKB:
		return bytesize / KB; // ת��ΪKB
	case SizeUnitMB:
		return bytesize / MB; // ת��ΪMB
	case SizeUnitGB:
		return bytesize / GB; // ת��ΪGB
	case SizeUnitTB:
		return bytesize / TB; // ת��ΪTB
	case SizeUnitPB:
		return bytesize / PB; // ת��ΪPB
	case SizeUnitEB:
		return bytesize / EB; // ת��ΪEB
	case SizeUnitZB:
		return bytesize / ZB; // ת��ΪZB
	case SizeUnitYB:
		return bytesize / YB; // ת��ΪYB
	default:
		return bytesize; // Ĭ�Ϸ����ֽ���
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
