#include "Vector.h"
#include "Matrix.h"

int main(int argc, char * argv[])
{
	EM::Matrix<float, 3, 3> A{9,-9,7,10,8,6,8,-9,9};
	EM::Matrix<float, 3, 3> B{1,-2,3,5,7,12,4,-19,-4};

	EM::Matrix<float, 3, 4> C{ 1,2,3,4,5,6,7,8,9,10,11,12 };
	EM::Matrix<float, 4, 3> D{ 1,2,3,4,5,6,7,8,9,10,11,12 };

	std::cout << C* D << std::endl;


	return 0;
}