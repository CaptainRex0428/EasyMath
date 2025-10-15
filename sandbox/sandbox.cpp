#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Vector<float,4> vectorA{1,2,3,4};
	std::cout << vectorA << std::endl;

	vectorA.x = 5;
	vectorA.zy = EM::Vector<float,2> {6,7};
	std::cout << vectorA << std::endl;

	std::cout << vectorA.xyz << std::endl;
	
	return 0;
}