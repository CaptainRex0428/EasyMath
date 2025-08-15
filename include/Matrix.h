#pragma once

#include "EasyMathAPI.h"
#include "Common.h"
#include "Vector.h"

namespace EM
{
	template <typename T, size_t rows, size_t cols>
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
			assert(InitializeList.size() == rows * cols && "Invalid initializer list size");
			std::copy(InitializeList.begin(), InitializeList.end(), data.begin());
		}

		// 访问元素（行主序）
		T& operator[](size_t idx)
		{
			limit_max(idx, rows * cols - 1);
			return data[idx];
		}

		const T& operator[](size_t idx) const
		{
			limit_max(idx, rows * cols - 1);
			return data[idx];
		}

		T& operator()(size_t row, size_t col)
		{
			limit_max(row, rows - 1);
			limit_max(col, cols - 1);

			return data[row * cols + col];
		}

		const T& operator()(size_t row, size_t col) const
		{
			limit_max(row, rows - 1);
			limit_max(col, cols - 1);

			return data[row * cols + col];
		}

		// 获取矩阵维度
		static constexpr size_t Rows() { return rows; }
		static constexpr size_t Cols() { return cols; }

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

	private:
		std::array<T, rows* cols> data;
	};

	typedef Matrix<float, 3, 3> Matrix3x3;
	typedef Matrix<float, 4, 4> Matrix4x4;
	typedef Matrix<float, 1, 4> Matrix1x4;
	typedef Matrix<float, 4, 1> Matrix4x1;

	template<typename T, size_t R, size_t C>
	Matrix<T, R, C> operator+(
		const Matrix<T, R, C>& a,
		const Matrix<T, R, C>& b)
	{
		Matrix<T, R, C> result;
		for (size_t i = 0; i < R * C; ++i) {
			result.data[i] = a.data[i] + b.data[i];
		}
		return result;
	}

	template<typename T, size_t R1, size_t C1, size_t C2>
	Matrix<T, R1, C2> operator*(
		const Matrix<T, R1, C1>& a,
		const Matrix<T, C1, C2>& b)
	{
		Matrix<T, R1, C2> result;
		for (size_t i = 0; i < R1; ++i) {
			for (size_t k = 0; k < C2; ++k) {
				T sum = 0;
				for (size_t j = 0; j < C1; ++j) {
					sum += a(i, j) * b(j, k);
				}
				result(i, k) = sum;
			}
		}
		return result;
	}

	template<typename T, size_t R, size_t C>
	Matrix<T, R, C> operator*(T scalar, const Matrix<T, R, C>& mat) {
		Matrix<T, R, C> result;
		for (size_t i = 0; i < R * C; ++i)
		{
			result.data[i] = scalar * mat.data[i];
		}
		return result;
	}

	Matrix4x1 operator*(Matrix4x4& matrix, Vector3& vector)
	{
		Matrix4x1 vectorM{ vector[x],vector[y],vector[z],1 };
		Matrix4x1 result = matrix * vectorM;

		return result;
	};

	Matrix4x1 operator*(Matrix4x4& matrix, Vector4& vector)
	{

		Matrix4x1 vectorM{ vector[x],vector[y],vector[z],vector[w] };
		Matrix4x1 result = matrix * vectorM;

		return result;

	};

	// 单位矩阵
	template<typename T, size_t N>
	Matrix<T, N, N> MTXIdentity()
	{
		static_assert(N > 0, "Identity matrix size must be positive");

		Matrix<T, N, N> mat;
		for (size_t i = 0; i < N; ++i) {
			mat(i, i) = 1;
		}
		return mat;
	}


	// 矩阵转置
	template<typename T, size_t R, size_t C>
	Matrix<T, C, R> MTXTranspose(const Matrix<T, R, C>& mat)
	{
		Matrix<T, C, R> result;
		for (size_t i = 0; i < R; ++i) {
			for (size_t j = 0; j < C; ++j) {
				result(j, i) = mat(i, j);
			}
		}
		return result;
	}

	// 3D旋转矩阵
	template<typename T>
	Matrix<T, 4, 4> MTXRotationX(T radians) {
		T cos = std::cos(radians);
		T sin = std::sin(radians);
		return {
			1,   0,    0, 0,
			0, cos, -sin, 0,
			0, sin,  cos, 0,
			0,   0,    0, 1
		};
	}

	template<typename T>
	Matrix<T, 4, 4> MTXRotationY(T radians) {
		T cos = std::cos(radians);
		T sin = std::sin(radians);
		return {
			cos , 0, sin, 0,
			0   , 1,   0, 0,
			-sin, 0, cos, 0,
			0   , 0,   0, 1
		};
	}

	template<typename T>
	Matrix<T, 4, 4> MTXRotationZ(T radians) {
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
	template<typename T>
	Matrix<T, 4, 4> MTXTranslation(T x, T y, T z)
	{
		Matrix<T, 4, 4> mat = MTXIdentity<T, 4>();
		mat(0, 3) = x;
		mat(1, 3) = y;
		mat(2, 3) = z;
		return mat;
	}

	// 3D缩放矩阵
	template<typename T>
	Matrix<T, 4, 4> MTXScale(T x, T y, T z)
	{
		Matrix<T, 4, 4> mat = MTXIdentity<T, 4>();
		mat(0, 0) = x;
		mat(1, 1) = y;
		mat(2, 2) = z;
		return mat;
	}
}