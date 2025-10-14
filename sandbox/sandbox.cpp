#include "Vector.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Color.h"

#define PRINT(var) std::cout << var << std::endl;

int main(int argc, char * argv[])
{
	EM::Vector<float,4> vector{ 1, 2, 3, 5};
	
	PRINT(vector.xyw)
	
	return 0;
}