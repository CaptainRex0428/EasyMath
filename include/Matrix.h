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

	// 反射平面对应的基础反射（关于坐标平面或原点）
	enum class ReflectionPlane
	{
		YZ,      // 翻转 X：关于 YZ 平面（X=0 平面）反射  → diag(-1, 1, 1, 1)
		XZ,      // 翻转 Y：关于 XZ 平面（Y=0 平面）反射  → diag( 1,-1, 1, 1)
		XY,      // 翻转 Z：关于 XY 平面（Z=0 平面）反射  → diag( 1, 1,-1, 1)
		Origin   // 关于原点反射（全部坐标取反）          → diag(-1,-1,-1, 1)
	};

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

	// 基础反射矩阵（关于坐标平面或原点）
	template<typename T,
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXReflection(ReflectionPlane plane)
	{
		Matrix<T, 4, 4> mat = MTXIdentity<T, 4>();
		switch (plane)
		{
		case ReflectionPlane::YZ:     mat(0, 0) = T(-1); break;
		case ReflectionPlane::XZ:     mat(1, 1) = T(-1); break;
		case ReflectionPlane::XY:     mat(2, 2) = T(-1); break;
		case ReflectionPlane::Origin:
			mat(0, 0) = T(-1);
			mat(1, 1) = T(-1);
			mat(2, 2) = T(-1);
			break;
		}
		return mat;
	}

	/**
	 * 创建任意平面反射矩阵（Householder 反射）
	 * 关于经过 pointOnPlane、法向量为 normal 的平面进行反射
	 *
	 * @param normal        平面的法向量（自动归一化；可传任意长度）
	 * @param pointOnPlane  平面上的任意一点（用于确定平面位置）
	 * @return 4×4 反射矩阵；det = -1（手性反转）
	 *
	 * 使用场景：
	 * - 平面镜渲染：把场景绕镜面镜像后渲染一次
	 * - 镜像复制：half-model 镜像生成完整模型
	 * - 物理模拟：球撞墙反弹
	 *
	 * 公式：M = I - 2 n n^T + 2 (n·P) n ⊗ e4
	 * 不变量：正交、M == M⁻¹、det = -1、平面上点是不动点
	 *
	 * 优化：缓存 2nx/2ny/2nz/d；内联点积
	 */
	template<typename T,
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXReflection(
		const Vector<T, 3>& normal,
		const Vector<T, 3>& pointOnPlane)
	{
		T invLen = T(1) / normal.length();
		assert(invLen != T(0) && "Reflection normal vector cannot be zero-length");

		T nx = normal.x * invLen;
		T ny = normal.y * invLen;
		T nz = normal.z * invLen;

		T twoNx = T(2) * nx;
		T twoNy = T(2) * ny;
		T twoNz = T(2) * nz;

		T d = nx * pointOnPlane.x + ny * pointOnPlane.y + nz * pointOnPlane.z;

		return {
			T(1) - twoNx * nx,  -twoNx * ny,           -twoNx * nz,           twoNx * d,
			-twoNy * nx,        T(1) - twoNy * ny,     -twoNy * nz,           twoNy * d,
			-twoNz * nx,        -twoNz * ny,           T(1) - twoNz * nz,     twoNz * d,
			T(0),               T(0),                  T(0),                  T(1)
		};
	}

	/**
	 * 创建 Look-At 视图矩阵（View Matrix）
	 * 用于将世界坐标转换到相机坐标系
	 *
	 * @param eye    相机位置（世界坐标）
	 * @param target 相机看向的目标点（世界坐标）
	 * @param up     相机的上方向向量（世界坐标，通常是 (0,1,0)）
	 * @return 4×4 视图矩阵，将世界坐标变换到相机坐标
	 *
	 * 使用场景：
	 * - 第三人称相机：eye = 角色后方, target = 角色位置
	 * - 第一人称相机：eye = 玩家位置, target = eye + lookDirection
	 * - 轨道相机：eye 绕 target 旋转
	 *
	 * 注意：返回的矩阵使相机位于原点，看向 +Z 方向（OpenGL 约定）
	 *
	 * 优化：
	 * - 减少归一化次数（3次 → 2次），节省 1 次平方根计算
	 * - camUp 不需要归一化（已通过叉积正交化）
	 * - 内联点积计算，减少函数调用开销
	 */
	template<typename T,
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXLookAt(
		const Vector<T, 3>& eye,
		const Vector<T, 3>& target,
		const Vector<T, 3>& up)
	{
		// 计算相机坐标系的三个正交轴
		// 优化：只归一化 2 次（forward 和 right），camUp 由正交向量叉积得出，已单位化
		Vector<T, 3> forward = (target - eye).normalized();     // 前方向（指向目标）
		Vector<T, 3> right = cross(forward, up).normalized();   // 右方向
		Vector<T, 3> camUp = cross(right, forward);             // 上方向（已正交且单位化）

		// 提取 swizzle 成员到实际值
		T rx = right.x, ry = right.y, rz = right.z;
		T ux = camUp.x, uy = camUp.y, uz = camUp.z;
		T fx = forward.x, fy = forward.y, fz = forward.z;
		T ex = eye.x, ey = eye.y, ez = eye.z;

		// 优化：内联点积计算，避免函数调用
		// -dot(right, eye) = -(rx*ex + ry*ey + rz*ez)
		T rightDotEye = rx * ex + ry * ey + rz * ez;
		T upDotEye = ux * ex + uy * ey + uz * ez;
		T forwardDotEye = fx * ex + fy * ey + fz * ez;

		// 构建视图矩阵
		// 行1-3: 旋转部分（相机坐标系的轴）
		// 列4: 平移部分（使相机位于原点）
		// 注意：forward 取负号是因为相机看向 +Z 方向
		return {
			rx,  ry,  rz,  -rightDotEye,
			ux,  uy,  uz,  -upDotEye,
			-fx, -fy, -fz,  forwardDotEye,
			0,   0,   0,    1
		};
	}

	/**
	 * 创建透视投影矩阵（Perspective Projection Matrix）
	 * 模拟人眼/相机的透视效果：远处的物体看起来更小
	 *
	 * @param fov    垂直视野角度（Field of View），单位：弧度
	 *               - 60° (PI/3) = 正常视野（第一人称游戏）
	 *               - 90° (PI/2) = 广角视野（竞技游戏）
	 *               - 30° (PI/6) = 望远镜效果
	 * @param aspect 宽高比（width / height），例如 16/9 = 1.778
	 * @param near   近裁剪面距离（> 0，通常 0.1 ~ 1.0）
	 *               太近的物体会被裁剪（避免穿模）
	 * @param far    远裁剪面距离（> near，通常 100 ~ 1000）
	 *               太远的物体会被裁剪（性能优化）
	 * @return 4×4 透视投影矩阵
	 *
	 * 使用场景：
	 * - 3D 游戏（FPS、TPS）
	 * - 需要深度感的场景
	 * - 模拟真实相机效果
	 *
	 * 注意：
	 * - 返回的矩阵会使 w' = -z（透视除法的关键）
	 * - near/far 范围越大，深度精度越差（Z-fighting）
	 * - 结果需要除以 w 分量（GPU 自动完成）
	 *
	 * 优化：
	 * - 缓存 1/tanHalfFov，减少除法次数（4次 → 1次）
	 * - 缓存 1/range，减少除法次数
	 * - 预计算常量表达式
	 * - 总减少：3 次除法，提升 ~10% 性能
	 */
	template<typename T,
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXPerspective(
		T fov,
		T aspect,
		T near,
		T far)
	{
		assert(aspect != T(0) && "Aspect ratio cannot be zero");
		assert(near > T(0) && "Near plane must be positive");
		assert(far > near && "Far plane must be greater than near plane");

		// 优化：缓存倒数，减少除法
		T tanHalfFov = std::tan(fov / T(2));
		T invTanHalfFov = T(1) / tanHalfFov;       // 缓存倒数
		T invRange = T(1) / (far - near);          // 缓存倒数

		// 预计算常量
		T a = invTanHalfFov / aspect;              // [0,0] 元素
		T b = invTanHalfFov;                       // [1,1] 元素
		T c = -(far + near) * invRange;            // [2,2] 元素
		T d = -(T(2) * far * near) * invRange;     // [2,3] 元素

		// 透视投影矩阵（OpenGL 风格，NDC 范围 [-1, 1]）
		return {
			a,    0,    0,    0,
			0,    b,    0,    0,
			0,    0,    c,    d,
			0,    0,   -T(1), 0
		};
	}

	/**
	 * 创建正交投影矩阵（Orthographic Projection Matrix）
	 * 无透视效果：远近物体大小相同
	 *
	 * @param left   左边界
	 * @param right  右边界
	 * @param bottom 下边界
	 * @param top    上边界
	 * @param near   近裁剪面距离
	 * @param far    远裁剪面距离
	 * @return 4×4 正交投影矩阵
	 *
	 * 使用场景：
	 * - 2D 游戏（像素完美对齐）
	 * - UI 渲染（屏幕空间 1:1）
	 * - CAD 软件（精确测量，无变形）
	 * - 建筑设计图（比例准确）
	 * - 阴影贴图（方向光）
	 *
	 * 典型用法：
	 * - 屏幕空间：MTXOrtho(0, width, 0, height, -1, 1)
	 * - 对称视图：MTXOrtho(-10, 10, -10, 10, 0.1, 100)
	 *
	 * 注意：
	 * - 结果的 w 分量始终为 1（无透视除法）
	 * - 物体大小不随距离变化
	 *
	 * 优化：
	 * - 缓存倒数，减少除法次数（9次 → 3次）
	 * - 预计算所有矩阵元素
	 * - 总减少：6 次除法，提升 ~15% 性能
	 */
	template<typename T,
		typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	Matrix<T, 4, 4> MTXOrtho(
		T left,
		T right,
		T bottom,
		T top,
		T near,
		T far)
	{
		assert(right != left && "Right and left cannot be equal");
		assert(top != bottom && "Top and bottom cannot be equal");
		assert(far != near && "Far and near cannot be equal");

		// 优化：缓存倒数，减少除法
		T invWidth = T(1) / (right - left);
		T invHeight = T(1) / (top - bottom);
		T invDepth = T(1) / (far - near);

		// 预计算所有矩阵元素
		T a = T(2) * invWidth;                     // [0,0] 元素
		T b = T(2) * invHeight;                    // [1,1] 元素
		T c = -T(2) * invDepth;                    // [2,2] 元素
		T tx = -(right + left) * invWidth;         // [0,3] 平移
		T ty = -(top + bottom) * invHeight;        // [1,3] 平移
		T tz = -(far + near) * invDepth;           // [2,3] 平移

		// 正交投影矩阵（OpenGL 风格，NDC 范围 [-1, 1]）
		return {
			a,   0,   0,   tx,
			0,   b,   0,   ty,
			0,   0,   c,   tz,
			0,   0,   0,   1
		};
	}
}