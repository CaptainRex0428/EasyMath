#include "Vector.h"
#include "Matrix.h"

int main(int argc, char * argv[])
{
	EM::Matrix<float, 4, 4> a{9,-9,7,10,8,6,8,-9,-7,2,-6,7,7,8,2,-8};

	std::cout << a << std::endl;

	std::cout << a.determinant() << std::endl;

	return 0;
}