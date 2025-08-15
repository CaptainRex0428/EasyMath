#pragma once

#ifdef EASYMATH_DLL
#define EASYMATH_API __declspec(dllexport)
#else
#ifdef EASYMATH_LIB
#define EASYMATH_API
#else
#define EASYMATH_API __declspec(dllimport)
#endif
#endif