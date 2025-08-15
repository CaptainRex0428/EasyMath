#pragma once

#include "EasyMathAPI.h"

#include <tuple>

namespace EM
{
	enum SystemSizeUnit
	{
		SizeUnitByte = 0,
		SizeUnitKB = 1,
		SizeUnitMB = 2,
		SizeUnitGB = 3,
		SizeUnitTB = 4,
		SizeUnitPB = 5,
		SizeUnitEB = 6,
		SizeUnitZB = 7,
		SizeUnitYB = 8
	};

	// transfer Byte size to specific unit output in double
	EASYMATH_API double ByteSizeTo(int bytesize, SystemSizeUnit unit);

	// transfer Byte size to a proper unit output in double
	EASYMATH_API std::tuple<const char*, double> ByteSizeConvert(uintmax_t bytesize);
}
