#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Vector<float, 3> vectorA{1,2,3};
	EM::Vector<float, 4> vectorB{4,5,6,7};

	std::cout << vectorA.xy << std::endl;
	return 0;
}