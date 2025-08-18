#pragma once

// ¼ì²âC++°æ±¾
#if __cplusplus >= 202002L
#define CPP20_OR_LATER
#elif __cplusplus >= 201703L
#define CPP17_OR_LATER
#elif __cplusplus >= 201402L
#define CPP14_OR_LATER
#elif __cplusplus >= 201103L
#define CPP11_OR_LATER
#else
#define CPP98_OR_EARLIER
#endif
