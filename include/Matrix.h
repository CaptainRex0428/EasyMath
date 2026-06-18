#pragma once

#include "EasyMathAPI.h"
#include "Common.h"
#include "Vector.h"

#include <iostream>
#include <array>
#include <initializer_list>
#include <cassert>
#include <iomanip>
#include <sstream>

namespace EasyMath
{
	template<typename T, size_t dimension, typename>
	class Vector;

	template <typename T, size_t rows, size_t cols, 
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	class Matrix
	{

	public:
		static_assert(rows > 0 && cols > 0, "Matrix row & column must be positive.");

		Matrix()
			:data({})
		{
		}

		Matrix(std::initializer_list<T> InitializeList)
		{
			assert(InitializeList.size() == rows * cols && "Invalid initializer list size.");

			std::copy(InitializeList.begin(), InitializeList.end(), data.begin());
		}

		Matrix(std::initializer_list<Vector<T,cols>> InitializeList)
		{
			assert(InitializeList.size() == rows && "Invalid initializer list size.");
    
			size_t row_idx = 0;
			for (const auto& vec : InitializeList)
			{
				for (size_t col = 0; col < cols; ++col)
				{
					data[row_idx * cols + col] = vec[col];
				}
				row_idx++;
			}
		}

		template<typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>
		Matrix(const Matrix<T2,rows,cols>& other):data({})
		{
			for (size_t i = 0; i < rows * cols; ++i)
			{
				data[i] = static_cast<T>(other[i]);
			}
		}

		virtual ~Matrix() = default;


/*
* 输出 
*/
		virtual std::ostream& print(std::ostream& out) const
		{
			out << std::fixed << std::setprecision(3);  // 设置小数精度为3，并使用固定点表示

			// 根据列宽决定每列的宽度
			size_t max_width = 0;
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < cols; ++j)
				{
					// 获取每个元素的字符串宽度
					std::ostringstream oss;
					oss << std::fixed << std::setprecision(3) << (*this)(i, j);  // 设置精度为3
					std::string str = oss.str();
					size_t element_width = str.length();
					max_width = std::max(max_width, element_width);  // 获取最大宽度
				}
			}

			std::cout << "Matrix " << rows << "x" << cols << std::endl;

			out << "-";
			out << std::setw(max_width * cols + cols+2);
			out << "-\n";

			for (size_t i = 0; i < rows; ++i)
			{
				out << "|";
				for (size_t j = 0; j < cols; ++j)
				{
					out << std::setw(max_width) << (*this)(i, j) << " ";
				}

				out << "|";
				out << "\n";  // 打印每一行后换行
			}
			
			out << "-";
			out << std::setw(max_width * cols + cols+2);
			out << "-\n";

			return out;
		
		}

/*
* 访问元素（行主序）
*/
		
		virtual T& operator[](size_t idx)
		{
			assert(idx < rows * cols && idx >= 0 && "Element index out of bound exception.");

			return data[idx];
		}

		virtual const T& operator[](size_t idx) const
		{
			assert(idx < rows * cols && idx >= 0 && "Element index out of bound exception.");
			
			return data[idx];
		}

		virtual T& operator()(size_t row, size_t col)
		{
			assert(row < rows && row >= 0 && "Element index out of bound exception.");
			assert(col < cols && col >= 0 && "Element index out of bound exception.");

			return data[row * cols + col];
		}

		virtual const T& operator()(size_t row, size_t col) const
		{
			assert(row < rows && row >= 0 && "Element index out of bound exception.");
			assert(col < cols && col >= 0 && "Element index out of bound exception.");

			return data[row * cols + col];
		}

		virtual T* Data() noexcept { return data.data(); }
		virtual const T* Data() const noexcept { return data.data(); }

		virtual T* begin() noexcept { return data.data(); }
		virtual const T* begin() const noexcept { return data.data(); }

		virtual T* end() noexcept { return data.data() + Size(); }
		virtual const T* end() const noexcept { return data.data() + Size(); }
		
		static constexpr size_t Size() noexcept { return rows*cols; }

		virtual T& at(size_t idx)
		{
			if (idx >= rows*cols) {
				throw std::out_of_range("Vector index out of range");
			}
			return data[idx];
		}

		virtual const T& at(size_t idx) const
		{
			if (idx >= rows*cols) {
				throw std::out_of_range("Vector index out of range");
			}
			return data[idx];
		}

		// 获取矩阵维度
		static constexpr size_t Rows() { return rows; }
		static constexpr size_t Cols() { return cols; }


/*
* 成员函数
*/
		// 创建去掉指定行列的子矩阵
		Matrix<T, rows - 1, cols - 1> submatrix(size_t row, size_t col) const
		{
			Matrix<T, rows - 1, cols - 1> submat;
			size_t sub_i = 0;
			for (size_t i = 0; i < rows; ++i)
			{
				if (i == row) continue;
				size_t sub_j = 0;
				for (size_t j = 0; j < cols; ++j)
				{
					if (j == col) continue;
					submat(sub_i, sub_j) = data[i * cols + j];
					++sub_j;
				}
				++sub_i;
			}
			return submat;
		}

		// 计算矩阵转置
		Matrix<T, cols, rows> transpose() const
		{
			Matrix<T, cols, rows> result{};

			for (size_t i = 0; i < rows; ++i) {
				for (size_t j = 0; j < cols; ++j) {
					result(j, i) = (*this)(i, j);
				}
			}
			return result;
		}

		
		// 计算行列式
		template<size_t R = rows, size_t C = cols>
		typename std::enable_if_t<R == C, T> determinant() const 
		{
			if constexpr (rows == 1) {
				// 1x1矩阵的行列式就是它的唯一元素
				return data[0];
			}
			else if constexpr (rows == 2)
			{
				// 2x2矩阵的行列式公式：ad - bc
				return data[0] * data[3] - data[1] * data[2];
			}
			else {
				T det = 0;
				for (size_t i = 0; i < cols; ++i) {
					// 根据行列式展开公式计算
					// 使用 submatrix() 函数生成子矩阵
					auto submat = submatrix(0, i);
					det += (i % 2 == 0 ? 1 : -1) * data[i] * submat.determinant();
				}
				return det;
			}
		}

		
		// 计算指定位置的代数余子式
		template<size_t R = rows, size_t C = cols>
		typename std::enable_if_t<R == C && R >= 2, T> cofactor(size_t row, size_t col) const
		{
			// 代数余子式 = (-1)^(i+j) * det(M_ij)
			// 其中 M_ij 是去掉第i行第j列后的子矩阵的行列式
			auto submat = submatrix(row, col);
			T minor = submat.determinant();

			// 计算符号：(-1)^(row+col)
			T sign = T(((row + col) % 2 == 0) ? 1 : -1);

			return sign * minor;
		}

		// 计算代数余子式矩阵（伴随矩阵的转置）
		template<size_t R = rows, size_t C = cols>
		typename std::enable_if_t<R == C, Matrix<T, rows, cols>> cofactorMatrix() const
		{
			static_assert(rows == cols, "Matrix must be square to compute cofactor matrix.");

			Matrix<T, rows, cols> result;

			for (size_t i = 0; i < rows; ++i) {
				for (size_t j = 0; j < cols; ++j) {
					result(i, j) = cofactor(i, j);
				}
			}

			return result;
		}

		// 计算伴随矩阵（代数余子式矩阵的转置）
		template<size_t R = rows, size_t C = cols>
		typename std::enable_if_t<R == C, Matrix<T, rows, cols>> adjugate() const
		{
			static_assert(rows == cols, "Matrix must be square to compute adjugate matrix.");

			auto cofMat = cofactorMatrix();
			return cofMat.transpose();
		}
		
		// 计算逆矩阵（使用伴随矩阵方法）
		template<size_t R = rows, size_t C = cols>
		typename std::enable_if_t<R == C, Matrix<T, rows, cols>> inverse() const
		{
			static_assert(rows == cols, "Matrix must be square to compute inverse.");

			T det = determinant();
			assert(det != 0 && "Matrix is singular (determinant is zero), cannot compute inverse.");

			auto adj = adjugate();
			return  (1.0 / det) * adj;
		}


/*
* 成员函数：复合赋值运算符
*
*/
		
		virtual Matrix<T, rows, cols>& operator+=(const Matrix<T, rows, cols>& matrix)
		{
			for (size_t i = 0; i < rows * cols; ++i)
			{
				data[i] += matrix.data[i];
			}
			return *this;
		}

		virtual Matrix<T, rows, cols>& operator+=(const T& scalar)
		{
			for (size_t i = 0; i < rows * cols; ++i)
			{
				data[i] += scalar;
			}
			return *this;
		}
		
		virtual Matrix<T, rows, cols>& operator-=(const Matrix<T, rows, cols>& matrix)
		{
			for (size_t i = 0; i < rows * cols; ++i)
			{
				data[i] -= matrix.data[i];
			}
			return *this;
		}

		virtual Matrix<T, rows, cols>& operator-=(const T& scalar)
		{
			for (size_t i = 0; i < rows * cols; ++i)
			{
				data[i] -= scalar;
			}
			return *this;
		}

		virtual Matrix<T, rows, cols>& operator*=(const Matrix<T, cols, cols>& matrix)
		{
			*this = *this * matrix;
			return *this;
		}
		
		virtual Matrix<T, rows, cols>& operator*=(const T& scalar)
		{
			for (size_t i = 0; i < rows * cols; ++i)
			{
				data[i] *= scalar;
			}
			return *this;
		}

		virtual Matrix<T, rows, cols>& operator/=(const T& scalar)
		{
			for (size_t i = 0; i < rows * cols; ++i)
			{
				data[i] /= scalar;
			}
			return *this;
		}
	private:
		

		friend Matrix<T, rows, cols> operator*(const Matrix<T, rows, cols>& matrix, const T& scalar)
		{
			Matrix<T, rows, cols> result;
			
			for (size_t i = 0; i < rows * cols; ++i) {
				result.data[i] = matrix.data[i] * static_cast<T>(scalar);
			}

			return result;
		}

		friend Matrix<T, rows, cols> operator*(const T& scalar, const Matrix<T, rows, cols>& matrix)
		{
			return matrix * scalar;
		}

		// 矩阵乘法运算符
		template<size_t cols2>
		friend Matrix<T, rows, cols2> operator*(const Matrix<T, rows, cols>& lhs, const Matrix<T, cols, cols2>& rhs)
		{
			Matrix<T, rows, cols2> result;

			for (size_t i = 0; i < rows; ++i) 
			{  
				for (size_t j = 0; j < cols2; ++j) 
				{
					T sum = T{ 0 };
					for (size_t k = 0; k < cols; ++k) 
					{
						sum += lhs(i, k) * rhs(k, j);
					}
					result(i, j) = sum;
				}
			}

			return result;
		}

		// 向量矩阵乘法运算
		// 矩阵 × 向量（列向量）: Matrix<T, rows, cols> × Vector<T, cols> -> Vector<T, rows>

		friend Matrix<T, rows, 1> operator*(const Matrix<T, rows, cols>&  matrix, const Vector<T, cols>& vec)
		{
			Matrix<T, rows, 1> result;

			for (size_t i = 0; i < rows; ++i)
			{
				T sum = T{ 0 };
				for (size_t j = 0; j < cols; ++j)
				{
					sum += matrix(i, j) * vec[j];
				}
				result[i] = sum;
			}

			return result;
		}

		// 向量（行向量）× 矩阵: Vector<T, rows> × Matrix<T, rows, cols> -> Vector<T, cols>

		friend Matrix<T, 1, cols> operator*(const Vector<T, cols>& vec, const Matrix<T, rows, cols>& matrix)
		{
			Matrix<T, 1, cols> result;

			for (size_t j = 0; j < cols; ++j)
			{
				T sum = T{ 0 };
				for (size_t i = 0; i < rows; ++i)
				{
					sum += vec[i] * matrix(i, j);
				}
				result[j] = sum;
			}

			return result;
		}

		// 矩阵加法运算符
		friend Matrix<T, rows, cols> operator+(const Matrix<T, rows, cols>& lhs, const Matrix<T, rows, cols>& rhs)
		{
			Matrix<T, rows, cols> result;

			// 逐元素相加
			for (size_t i = 0; i < rows * cols; ++i) {
				result.data[i] = lhs.data[i] + rhs.data[i];
			}

			return result;
		}

		friend Matrix<T, rows, cols> operator+(const Matrix<T, rows, cols>& lhs, const T& scalar)
		{
			Matrix<T, rows, cols> result;

			// 逐元素相加
			for (size_t i = 0; i < rows * cols; ++i) {
				result.data[i] = lhs.data[i] + scalar;
			}

			return result;
		}

		// 矩阵减法运算符
		friend Matrix<T, rows, cols> operator-(const Matrix<T, rows, cols>& lhs, const Matrix<T, rows, cols>& rhs)
		{
			Matrix<T, rows, cols> result;

			// 逐元素相减
			for (size_t i = 0; i < rows * cols; ++i) {
				result.data[i] = lhs.data[i] - rhs.data[i];
			}

			return result;
		}

		friend Matrix<T, rows, cols> operator-(const Matrix<T, rows, cols>& lhs, const T& scalar)
		{
			Matrix<T, rows, cols> result;

			// 逐元素相加
			for (size_t i = 0; i < rows * cols; ++i) {
				result.data[i] = lhs.data[i] - scalar;
			}

			return result;
		}


		// 矩阵取负运算符（一元减号）
		friend Matrix<T, rows, cols> operator-(const Matrix<T, rows, cols>& matrix)
		{
			Matrix<T, rows, cols> result;

			// 逐元素取负
			for (size_t i = 0; i < rows * cols; ++i) 
			{
				result.data[i] = -matrix.data[i];
			}

			return result;
		}

		template<typename T2, size_t rows2, size_t cols2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>
		friend bool operator==(const Matrix<T,rows,cols>& matrixA, const Matrix<T2,rows2,cols2>& matrixB)
		{
			if (matrixA.Rows()!= matrixB.Rows() || matrixA.Cols() != matrixB.Cols())
			{
				return false;
			}

			if constexpr (std::is_floating_point_v<T>)
			{
				for (size_t i = 0; i < matrixA.Rows(); ++i)
				{
					for (size_t j = 0; j < matrixA.Cols(); ++j)
					{
						if (std::abs(matrixA(i, j) - matrixB(i, j)) > T{ NEARZERO_THRESHOLD })
						{
							return false;
						}
					}
				}
			}
			else
			{
				for (size_t i = 0; i < matrixA.Rows(); ++i)
				{
					for (size_t j = 0; j < matrixA.Cols(); ++j)
					{
						if (matrixA(i, j) != matrixB(i, j))
						{
							return false;
						}
					}
				}
			}

			return true;
			
		}
		
		template<typename T2, size_t rows2, size_t cols2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>
		friend bool operator!=(const Matrix<T,rows,cols>& matrixA, const Matrix<T2,rows2,cols2>& matrixB)
		{
			return !(matrixA == matrixB);
		}

	private:
		std::array<T, rows* cols> data;

	};

	template<typename T, size_t rows, size_t cols>
		std::ostream& operator<<(std::ostream& out, const Matrix<T,rows, cols>& matrix)
	{
		return matrix.print(out);
	}
	
	// 单位矩阵
	template<typename T, size_t N, 
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, N, N> MTXIdentity()
	{
		static_assert(N > 0, "Identity matrix size must be positive");

		Matrix<T, N, N> mat;
		for (size_t i = 0; i < N; ++i) {
			mat(i, i) = 1;
		}
		return mat;
	}

	// 3D旋转矩阵
	template<typename T, 
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXRotationX(const T& radians) {
		T cos = std::cos(radians);
		T sin = std::sin(radians);
		return {
			1,   0,    0, 0,
			0, cos, -sin, 0,
			0, sin,  cos, 0,
			0,   0,    0, 1
		};
	}

	template<typename T, 
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXRotationY(const T& radians) {
		T cos = std::cos(radians);
		T sin = std::sin(radians);
		return {
			cos , 0, sin, 0,
			0   , 1,   0, 0,
			-sin, 0, cos, 0,
			0   , 0,   0, 1
		};
	}

	template<typename T,
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXRotationZ(const T& radians) {
		T cos = std::cos(radians);
		T sin = std::sin(radians);
		return {
			cos, -sin, 0, 0,
			sin, cos , 0, 0,
			0,   0,    1, 0,
			0,   0,    0, 1
		};
	}

	// 3D平移矩阵
	template<typename T,
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXTranslation(const T& x, const T& y, const T& z, bool usedWithOrient = false)
	{
		Matrix<T, 4, 4> mat = MTXIdentity<T, 4>();
		mat(0, 3) = x;
		mat(1, 3) = y;
		mat(2, 3) = z;

		if (usedWithOrient)
		{
			mat(3, 3) = 0;
		}

		return mat;
	}

	// 3D缩放矩阵
	template<typename T,
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXScale(const T& x, const T& y, const T& z)
	{
		Matrix<T, 4, 4> mat = MTXIdentity<T, 4>();
		mat(0, 0) = x;
		mat(1, 1) = y;
		mat(2, 2) = z;
		return mat;
	}
}