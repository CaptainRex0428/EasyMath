#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Vector<float,2> vectorA{ 1, 2};
	EM::Vector<float,3> vectorB{ 3, 4, 5};
	EM::Vector<float,4> vectorC{ 6, 7, 8, 9};
	
	std::cout << vectorB.toColMatrix() << std::endl;
	std::cout << vectorB.toRowMatrix() << std::endl;
	std::cout << vectorB.length() << std::endl;
	std::cout << vectorB.normalized() << std::endl;
	std::cout << vectorB.isNormalized() << std::endl;
	std::cout << vectorB.isZero() << std::endl;
	std::cout << vectorB.lerp(EM::Vector<float,3>{1,1,1},0.5f) << std::endl;
	std::cout << vectorB.reflect(EM::Vector<float,3>{1,1,1}.normalize()) << std::endl;
	std::cout << vectorB.project(EM::Vector<float,3>{1,1,1}.normalize()) << std::endl;
	std::cout << vectorB.toHomogeneous(1) << std::endl;
	std::cout << vectorB.toTranslationMatrix() << std::endl;
	std::cout << vectorB.skewSymmetric() << std::endl;
	
	std::cout << vectorA.skewSymmetric_2D() << std::endl;
	std::cout << vectorC.fromHomogeneous() << std::endl;
	
	
	
	return 0;
}