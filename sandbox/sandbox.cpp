#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Vector<float,3> vector{ 1, 2, 3};
	
	PRINT(vector.xz)
	
	return 0;
}