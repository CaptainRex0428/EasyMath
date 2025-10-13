#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Vector<float, 3> A{ 1, 2, 3};
	EM::Vector<float, 3> B{ 1, 2, 3};
	EM::Vector<float, 4> C{ 1, 2, 3, 5};

	PRINT(C)

	EM::Matrix<float, 3, 3> CM{ 1,3,3,-9,7,6,7,8,9 };
	EM::Matrix<float, 2, 2> DM{ 1,2,3,4};
	EM::Matrix<float, 4, 2> EM{ 1,2,3,4,5,6,7,8};

	PRINT(CM.inverse());
	PRINT(CM.transpose());
	PRINT(A.skewSymmetric());
	PRINT(C.skewSymmetric());
	PRINT(A.toTranslationMatrix());

	EM::Quaternion<float> QA(A,4);
	PRINT(QA.normalized())
	PRINT(QA)

	EM::Color<float> CA(.4f);
	PRINT(CA.luminance())
	
	return 0;
}