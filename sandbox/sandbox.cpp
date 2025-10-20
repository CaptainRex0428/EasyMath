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
	
	return 0;
}