#include "Vector.h"
#include "Matrix.h"

int main(int argc, char * argv[])
{
	EM::Vector<float, 5> A{ 1, 2, 3, 5, 6};
	EM::Vector<float, 5> B{ 1, 2, 3, 5, 6};

	std::cout << A.toRowMatrix() << std::endl;

	return 0;
}