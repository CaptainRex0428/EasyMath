#include <iostream>
#include <ostream>

#ifdef ENGINE_API
#include "EasyMath/include/Color.h"
#include "EasyMath/include/Vector.h"
#include "EasyMath/include/Matrix.h"
#include "EasyMath/include/Quaternion.h"
#else
#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"
#endif

#ifdef ENGINE_API

#elif
#define PRINT(var) std::cout << var << std::endl;
int main(int argc, char * argv[])
{

	EM::sRGBColor<float> ColorS{0.4f,0.6f,0.9f,1.f};
	EM::LinearColor<float> ColorL{0.4f,0.6f,0.9f,1.f};
	std::cout << ColorS.luminance() << std::endl;
	std::cout << ColorL.luminance()<< std::endl;

	EM::LinearColor<float> ColorL2 = ColorS.toLinear();
	EM::sRGBColor<float> ColorS2 = ColorL.toSRGB();
	
	std::cout << EM::sRGBColor<float>::fromHSV(EM::Vector<float,3>{1,1,1},1) << std::endl;
	std::cout << EM::sRGBColor<float>::fromHSL(EM::Vector<float,3>{1,1,1},1)<< std::endl;
	std::cout << EM::sRGBColor<float>::fromHSI(EM::Vector<float,3>{1,1,1},1)<< std::endl;
	std::cout << EM::LinearColor<float>::fromHSV(EM::Vector<float,3>{1,1,1},1)<< std::endl;
	std::cout << EM::LinearColor<float>::fromHSL(EM::Vector<float,3>{1,1,1},1)<< std::endl;
	std::cout << EM::LinearColor<float>::fromHSI(EM::Vector<float,3>{1,1,1},1)<< std::endl;
	
	std::cout << ColorL2 * ColorS2 << std::endl;

	
	
	return 0;
}

#endif
