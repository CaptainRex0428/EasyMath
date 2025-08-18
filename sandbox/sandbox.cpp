#include "Vector.h"
#include "Matrix.h"

int main(int argc, char * argv[])
{
	EM::Vector<float, 3> A{ 1, 2, 3};
	EM::Vector<float, 3> B{ 1, 2, 3};

	EM::Matrix<float, 3,3> CM{ 1,2,3,4,5,6,7,8,9};

	std::cout << B * CM << std::endl;
	std::cout << CM * B << std::endl;

	return 0;
}