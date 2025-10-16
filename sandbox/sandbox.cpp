#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Matrix<double,4,3> matrixA {1,2,3,4,5,6,7,8,9,10,11,12};
	EM::Matrix<double,4,4> matrixB {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	
	matrixA.submatrix(0,1);
	matrixA.transpose();

	matrixB.submatrix(1,2);
	matrixB.transpose();
	matrixB.determinant();
	matrixB.cofactorMatrix();
	matrixB.adjugate();
	matrixB.inverse();
	
	return 0;
}