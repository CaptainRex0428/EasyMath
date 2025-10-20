#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Vector<float, 3> vectorA{1,2,3};
	
	std::cout << vectorA.z + 5.f << std::endl;
	std::cout << vectorA.zy + 5.f << std::endl;
	std::cout << 5.f + vectorA.z  << std::endl;
	std::cout << 5.f + vectorA.zy  << std::endl;
	std::cout << vectorA.z - 5.f << std::endl;
	std::cout << vectorA.zy - 5.f << std::endl;
	std::cout << 5.f-vectorA.z << std::endl;
	//std::cout << 5.f-vectorA.zy << std::endl;
	std::cout << vectorA.z * 5.f << std::endl;
	std::cout << vectorA.zy * 5.f << std::endl;
	std::cout << 5.f * vectorA.z << std::endl;
	std::cout << 5.f * vectorA.zy << std::endl;
	std::cout << vectorA.z / 5.f << std::endl;
	std::cout << vectorA.zy / 5.f << std::endl;
	std::cout << 5.f / vectorA.z << std::endl;

	EM::Vector<float, 4> vectorB{4,5,6,7};
	std::cout << vectorB.w << std::endl;
	std::cout << vectorB.y * vectorB << std::endl;
	std::cout << vectorB / vectorB.w << std::endl;

	std::cout << vectorB.xyw * vectorA << std::endl;
	std::cout << vectorB.xyw / vectorA << std::endl;

	std::cout << vectorA * vectorB.xyw  << std::endl;
	std::cout << vectorA / vectorB.xyw << std::endl;

	std::cout << vectorA.x + vectorB.xyw << std::endl;
	std::cout << vectorA.x * vectorB.xyw << std::endl;
	std::cout << vectorB.xyw + vectorA.x << std::endl;
	std::cout << vectorB.xyw * vectorA.x << std::endl;
	std::cout << vectorB.xyw - vectorA.x << std::endl;
	std::cout << vectorB.xyw / vectorA.y << std::endl;

	std::cout << vectorB.xyw + vectorA.xyz << std::endl;
	std::cout << vectorB.xyw * vectorA.xyz << std::endl;
	std::cout << vectorB.xyw - vectorA.xyz << std::endl;
	std::cout << vectorB.xyw / vectorA.xyz << std::endl;

	std::cout << vectorB.xyw + vectorA.xyz << std::endl;
	std::cout << vectorB.xyw * vectorA.xyz << std::endl;
	std::cout << vectorB.xyw - vectorA.xyz << std::endl;
	std::cout << vectorB.xyw / vectorA.xyz << std::endl;

	std::cout << dot(vectorB.x, vectorA.y) << std::endl;
	std::cout << cross(vectorB.xyw, vectorA.xyz) << std::endl;
	std::cout << cross(vectorB.xy, vectorA.xy) << std::endl;

	std::cout << vectorB.x * vectorA.x << std::endl;

	std::cout << distance(vectorB.xy, vectorA.xy) << std::endl;

	EM::Color<float> cA{0.3f,0.5f,0.7f,1.f};
	EM::Color<float> cB{0.2f,0.8f,0.9f,1.f};

	std::cout << cA.luminance(EM::LuminanceStandard::Rec601) << std::endl;
	std::cout << cA.luminance(EM::LuminanceStandard::Rec709) << std::endl;
	std::cout << cB.luminance(EM::LuminanceStandard::Rec601) << std::endl;
	std::cout << cB.luminance(EM::LuminanceStandard::Rec709) << std::endl;

	cA.r = 1000.f;

	std::cout << cA.g << std::endl;


	EM::sRGBColor<float>::fromHSV(EM::Vector<float,3>{1,1,1},1);
	EM::sRGBColor<float>::fromHSL(EM::Vector<float,3>{1,1,1},1);
	EM::sRGBColor<float>::fromHSI(EM::Vector<float,3>{1,1,1},1);
	EM::LinearColor<float>::fromHSV(EM::Vector<float,3>{1,1,1},1);
	EM::LinearColor<float>::fromHSL(EM::Vector<float,3>{1,1,1},1);
	EM::LinearColor<float>::fromHSI(EM::Vector<float,3>{1,1,1},1);

	// EM::sRGBColor<float> sRGBA;
	// sRGBA.toLinear<>();
	//
	// EM::LinearColor<float> LinearA;
	// LinearA.toSRGB<>();
	//
	// EM::HSV<float> HSVA;
	// HSVA.toSRGB<>();
	// HSVA.toLinear<>();
	
	return 0;
}