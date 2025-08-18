#pragma once

#include "Vector.h"
#include "Matrix.h"

namespace EM
{
	// ========== Vector 到 Matrix 的转换 ==========

	// 将 Vector 转换为列向量矩阵
	template<typename T, size_t dimension>
	Matrix<T, dimension, 1> VectorCMatrix(const Vector<T, dimension>& vec)
	{
		Matrix<T, dimension, 1> result;
		for (size_t i = 0; i < dimension; ++i) {
			result(i, 0) = vec[i];
		}
		return result;
	}

	// 将 Vector 转换为行向量矩阵
	template<typename T, size_t dimension>
	Matrix<T, 1, dimension> VectorRMatrix(const Vector<T, dimension>& vec)
	{
		Matrix<T, 1, dimension> result;
		for (size_t i = 0; i < dimension; ++i) {
			result(0, i) = vec[i];
		}
		return result;
	}

	// 将多个 Vector 组合成矩阵（每个 Vector 作为一列）
	template<typename T, size_t dimension, size_t numVectors>
	Matrix<T, dimension, numVectors> VectorsToMatrix(const std::array<Vector<T, dimension>, numVectors>& vectors)
	{
		Matrix<T, dimension, numVectors> result;
		for (size_t col = 0; col < numVectors; ++col) {
			for (size_t row = 0; row < dimension; ++row) {
				result(row, col) = vectors[col][row];
			}
		}
		return result;
	}

	// 将多个 Vector 组合成矩阵（每个 Vector 作为一行）
	template<typename T, size_t dimension, size_t numVectors>
	Matrix<T, numVectors, dimension> VectorsToMatrixAsRows(const std::array<Vector<T, dimension>, numVectors>& vectors)
	{
		Matrix<T, numVectors, dimension> result;
		for (size_t row = 0; row < numVectors; ++row) {
			for (size_t col = 0; col < dimension; ++col) {
				result(row, col) = vectors[row][col];
			}
		}
		return result;
	}



	// ========== Matrix 到 Vector 的转换 ==========

	// 从列向量矩阵提取 Vector
	template<typename T, size_t rows>
	Vector<T, rows> ColumnMatrixToVector(const Matrix<T, rows, 1>& matrix)
	{
		Vector<T, rows> result;
		for (size_t i = 0; i < rows; ++i) {
			result[i] = matrix(i, 0);
		}
		return result;
	}

	// 从行向量矩阵提取 Vector
	template<typename T, size_t cols>
	Vector<T, cols> RowMatrixToVector(const Matrix<T, 1, cols>& matrix)
	{
		Vector<T, cols> result;
		for (size_t i = 0; i < cols; ++i) {
			result[i] = matrix(0, i);
		}
		return result;
	}

	// 从矩阵的指定列提取 Vector
	template<typename T, size_t rows, size_t cols>
	Vector<T, rows> ExtractColumn(const Matrix<T, rows, cols>& matrix, size_t col)
	{
		assert(col < cols && "Column index out of bounds");
		Vector<T, rows> result;
		for (size_t i = 0; i < rows; ++i) {
			result[i] = matrix(i, col);
		}
		return result;
	}

	// 从矩阵的指定行提取 Vector
	template<typename T, size_t rows, size_t cols>
	Vector<T, cols> ExtractRow(const Matrix<T, rows, cols>& matrix, size_t row)
	{
		assert(row < rows && "Row index out of bounds");
		Vector<T, cols> result;
		for (size_t i = 0; i < cols; ++i) {
			result[i] = matrix(row, i);
		}
		return result;
	}

	// 从矩阵的对角线提取 Vector
	template<typename T, size_t N>
	Vector<T, N> ExtractDiagonal(const Matrix<T, N, N>& matrix)
	{
		Vector<T, N> result;
		for (size_t i = 0; i < N; ++i) {
			result[i] = matrix(i, i);
		}
		return result;
	}

	// ========== Matrix-Vector 运算 ==========

	// 矩阵乘以向量（将向量视为列向量）
	template<typename T, size_t rows, size_t cols>
	Vector<T, rows> operator*(const Matrix<T, rows, cols>& matrix, const Vector<T, cols>& vec)
	{
		Vector<T, rows> result;
		for (size_t i = 0; i < rows; ++i) {
			T sum = T{ 0 };
			for (size_t j = 0; j < cols; ++j) {
				sum += matrix(i, j) * vec[j];
			}
			result[i] = sum;
		}
		return result;
	}

	// 向量乘以矩阵（将向量视为行向量）
	template<typename T, size_t rows, size_t cols>
	Vector<T, cols> operator*(const Vector<T, rows>& vec, const Matrix<T, rows, cols>& matrix)
	{
		Vector<T, cols> result;
		for (size_t j = 0; j < cols; ++j) {
			T sum = T{ 0 };
			for (size_t i = 0; i < rows; ++i) {
				sum += vec[i] * matrix(i, j);
			}
			result[j] = sum;
		}
		return result;
	}

	// ========== 特殊转换函数 ==========

	// 创建反对称矩阵（用于3D向量的叉积运算）
	template<typename T>
	Matrix<T, 3, 3> VectorToSkewSymmetricMatrix(const Vector<T, 3>& vec)
	{
		return {
			 T{0}, -vec[z],  vec[y],
			 vec[z],  T{0}, -vec[x],
			-vec[y],  vec[x],  T{0}
		};
	}

	// 创建齐次坐标变换矩阵（从3D向量创建平移矩阵）
	template<typename T>
	Matrix<T, 4, 4> VectorToTranslationMatrix(const Vector<T, 3>& translation)
	{
		Matrix<T, 4, 4> result = MTXIdentity<T, 4>();
		result(0, 3) = translation[x];
		result(1, 3) = translation[y];
		result(2, 3) = translation[z];
		return result;
	}

	// 从变换矩阵提取平移向量
	template<typename T>
	Vector<T, 3> ExtractTranslationFromMatrix(const Matrix<T, 4, 4>& matrix)
	{
		return { matrix(0, 3), matrix(1, 3), matrix(2, 3) };
	}

	// 从变换矩阵提取缩放向量
	template<typename T>
	Vector<T, 3> ExtractScaleFromMatrix(const Matrix<T, 4, 4>& matrix)
	{
		Vector<T, 3> scale;
		// 计算每列的长度来获取缩放因子
		for (size_t i = 0; i < 3; ++i) {
			T length = T{ 0 };
			for (size_t j = 0; j < 3; ++j) {
				length += matrix(j, i) * matrix(j, i);
			}
			scale[i] = std::sqrt(length);
		}
		return scale;
	}

	// ========== 便利函数 ==========

	// 创建外积矩阵（u ⊗ v）
	template<typename T, size_t M, size_t N>
	Matrix<T, M, N> OuterProduct(const Vector<T, M>& u, const Vector<T, N>& v)
	{
		Matrix<T, M, N> result;
		for (size_t i = 0; i < M; ++i) {
			for (size_t j = 0; j < N; ++j) {
				result(i, j) = u[i] * v[j];
			}
		}
		return result;
	}

	// 使用向量构建旋转矩阵（Rodrigues公式）
	template<typename T>
	Matrix<T, 3, 3> VectorToRotationMatrix(const Vector<T, 3>& axis, T angle)
	{
		Vector<T, 3> normalizedAxis = axis.normalized();
		T cos_theta = std::cos(angle);
		T sin_theta = std::sin(angle);
		T one_minus_cos = T{ 1 } - cos_theta;

		T x = normalizedAxis[0];
		T y = normalizedAxis[1];
		T z = normalizedAxis[2];

		return {
			cos_theta + x * x * one_minus_cos,     x * y * one_minus_cos - z * sin_theta,  x * z * one_minus_cos + y * sin_theta,
			y * x * one_minus_cos + z * sin_theta,  cos_theta + y * y * one_minus_cos,     y * z * one_minus_cos - x * sin_theta,
			z * x * one_minus_cos - y * sin_theta,  z * y * one_minus_cos + x * sin_theta,   cos_theta + z * z * one_minus_cos
		};
	}

	// ========== 类型别名 ==========

	// 常用的矩阵-向量组合类型
	using Matrix2f = Matrix<float, 2, 2>;
	using Matrix3f = Matrix<float, 3, 3>;
	using Matrix4f = Matrix<float, 4, 4>;

	using Matrix2d = Matrix<double, 2, 2>;
	using Matrix3d = Matrix<double, 3, 3>;
	using Matrix4d = Matrix<double, 4, 4>;

	// 便利的转换函数别名
	template<typename T>
	using Vec3ToMat4 = Matrix<T, 4, 4>(*)(const Vector<T, 3>&);

	template<typename T>
	using Mat4ToVec3 = Vector<T, 3>(*)(const Matrix<T, 4, 4>&);
}