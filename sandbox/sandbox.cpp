#include "Vector.h"
#include "Matrix.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Vector<float, 3> A{ 1, 2, 3};
	EM::Vector<float, 3> B{ 1, 2, 3};

	EM::Matrix<float, 3, 3> CM{ 1,3,3,-9,7,6,7,8,9 };
	EM::Matrix<float, 2, 2> DM{ 1,2,3,4};
	EM::Matrix<float, 4, 2> EM{ 1,2,3,4,5,6,7,8};

	PRINT(CM.inverse());
	PRINT(CM.transpose());

	return 0;
}