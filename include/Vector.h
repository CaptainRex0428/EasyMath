#pragma once

#include "EasyMathAPI.h"
#include "Common.h"
#include "Version.h"

#include <cmath>
#include <iostream>
#include <array>
#include <cassert>
#include <algorithm>
#include <vector>

namespace EasyMath
{
	template <typename T, size_t rows, size_t cols, typename>
	class Matrix;

	template<typename T, typename>
	class Quaternion;

	template <typename T, size_t dimension, int... Indices>
	struct Swizzle;

	template<typename T, size_t dimension, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	class Vector
	{
	public:
		static_assert(std::is_arithmetic_v<T>, "Vector element type must be arithmetic");
		static_assert(dimension > 0, "vector dimension must be positive.");
		
		Vector()
			:data({})
		{
		}

		explicit Vector(const T& value)
		{
			data.fill(value);
		}

		Vector(std::initializer_list<T> InitializeList)
		{
			assert(InitializeList.size() == dimension && "[Vector] Invalid initializer list size");

			std::copy(InitializeList.begin(), InitializeList.end(), data.begin());
		}

		template<typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>
		Vector(const Vector<T2, dimension>& other):data({})
		{
			for (size_t i = 0; i < dimension; ++i)
			{
				data[i] = static_cast<T>(other[i]);
			}
		}

		Vector(const Vector&) = default;

		Vector(Vector&&) = default;

		virtual ~Vector() = default;

/*
* 数据访问
*
*/
		virtual T& operator[](size_t idx)
		{
			assert(idx < dimension && "Index out of bounds");
			return data[idx];
		}

		virtual const T& operator[](size_t idx) const
		{
			assert(idx < dimension && "Index out of bounds");
			return data[idx];
		}
		
		virtual T* Data() noexcept { return data.data(); }
		virtual const T* Data() const noexcept { return data.data(); }

		virtual T* begin() noexcept { return data.data(); }
		virtual const T* begin() const noexcept { return data.data(); }

		virtual T* end() noexcept { return data.data() + dimension; }
		virtual const T* end() const noexcept { return data.data() + dimension; }
		
		static constexpr size_t Size() noexcept { return dimension; }
		virtual constexpr size_t getDimension() const noexcept { return dimension; }
		static constexpr size_t Dimension() noexcept { return dimension; }

		virtual T& at(size_t idx)
		{
			if (idx >= dimension) {
				throw std::out_of_range("Vector index out of range");
			}
			return data[idx];
		}

		virtual const T& at(size_t idx) const
		{
			if (idx >= dimension) {
				throw std::out_of_range("Vector index out of range");
			}
			return data[idx];
		}
		
/*
* 直接被operator<<调用
*
*/
		virtual std::ostream& print(std::ostream& out) const
		{
			std::string mode = "";

			switch (dimension)
			{
			case 1: mode = ": x"; break;
			case 2: mode = ": xy"; break;
			case 3: mode = ": xyz"; break;
			case 4: mode = ": xyzw"; break;
			default: break;
			}
			
			out << "(";
			for (size_t i = 0; i < dimension; ++i) 
			{
				out << (*this).data[i];
				if (i < dimension - 1) {
					out << ", ";
				}
			}
			out << ") (V mode" << mode << ")";
			return out;
		}

/*
* 成员函数：矩阵形式
*
*/
		
		virtual Matrix<T, dimension, 1, std::enable_if_t<std::is_arithmetic_v<T>>> toColMatrix()
		{
			Matrix<T, dimension, 1, std::enable_if_t<std::is_arithmetic_v<T>>> result{};
			for (size_t i = 0; i < dimension; ++i) 
			{
				result(i, 0) = this->data[i];
			}
			return result;
		}

		virtual Matrix<T, 1, dimension, std::enable_if_t<std::is_arithmetic_v<T>>> toRowMatrix()
		{
			Matrix<T, 1, dimension, std::enable_if_t<std::is_arithmetic_v<T>>> result{};
			for (size_t i = 0; i < dimension; ++i) 
			{
				result(0, i) = this->data[i];
			}
			return result;
		}

/*
* 成员函数：向量运算
*
*/
		virtual constexpr T lengthSquared(bool dimensionalityReduction = false) const noexcept
		{
			if (!dimensionalityReduction) {
				T sum = T(0);
				for (size_t idx = 0; idx < dimension; ++idx)
				{
					sum += data[idx] * data[idx];
				}
				return sum;
			}

			T sum = T(0);
			size_t maxDim = std::min(dimension, size_t(3));
			for (size_t idx = 0; idx < maxDim; ++idx)
			{
				sum += data[idx] * data[idx];
			}
			return sum;
		}

		virtual T length(bool dimensionalityReduction = false) const noexcept
		{
			return std::sqrt(lengthSquared(dimensionalityReduction));
		}

		virtual Vector<T, dimension> normalized(bool dimensionalityReduction = false) const
		{
			T len = length(dimensionalityReduction);

			if (len == T{ 0 })
			{
				return Vector<T, dimension>{};
			}

			Vector<T, dimension> result;

			if (dimensionalityReduction)
				{
				size_t maxDim = std::min(dimension, size_t(3));
				for (size_t idx = 0; idx < maxDim; ++idx)
				{
					result[idx] = data[idx] / len;
				}
				// 保持其他维度不变
				for (size_t idx = maxDim; idx < dimension; ++idx)
				{
					result[idx] = data[idx];
				}
			}
			else
			{
				for (size_t idx = 0; idx < dimension; ++idx)
				{
					result[idx] = data[idx] / len;
				}
			}

			return result;
		}
		
		virtual Vector<T, dimension> normalize(bool dimensionalityReduction = false)
		{
			*this = normalized(dimensionalityReduction);
			return *this;
		}

		virtual bool isNormalized(bool dimensionalityReduction = true, T epsilon = T{ NEARZERO_THRESHOLD }) const noexcept
		{
			return std::abs(length(dimensionalityReduction) - T{ 1 }) <= epsilon;
		}
		
		virtual bool isZero(T epsilon = T{ NEARZERO_THRESHOLD }) const noexcept
		{
			for (size_t i = 0; i < dimension; ++i) 
			{
				if (std::abs(data[i]) > epsilon) 
				{
					return false;
				}
			}
			return true;
		}
		
		// 线性插值
		virtual Vector<T, dimension> lerp(const Vector<T, dimension>& other, T t) const noexcept
		{
			Vector<T, dimension> result;
			for (size_t i = 0; i < dimension; ++i) {
				result.data[i] = data[i] + t * (other.data[i] - data[i]);
			}
			return result;
		}

		// 反射（需要基础单位法向量）
		virtual Vector<T, dimension> reflect(const Vector<T, dimension>& normal) const
		{
			T dot_product = dot(*this, normal.normalized());
			return T{ 2 } *dot_product * normal - * this;
		}

		// 投影
		virtual Vector project(const Vector& onto) const
		{
			T dot_product = dot(*this, onto);
			T onto_length_sq = onto.lengthSquared();
			if (onto_length_sq == T{ 0 }) 
			{
				return Vector{};
			}
			return (dot_product / onto_length_sq) * onto;
		}

		virtual Vector project(const Vector& onto, bool dimensionalityReduction) const
		{
			T dot_product = dot(*this, onto);
			T onto_length_sq = onto.lengthSquared(dimensionalityReduction);
			if (onto_length_sq == T{ 0 }) {
				return Vector{};
			}
			return (dot_product / onto_length_sq) * onto;
		}

		// 向量的齐次坐标相关操作

		template<size_t newDim = dimension + 1>
		Vector<T, newDim> toHomogeneous(T w = T{ 1 }) const
		{
			static_assert(newDim > dimension, "New dimension must be larger");
			Vector<T, newDim> result;

			for (size_t i = 0; i < dimension; ++i) {
				result[i] = data[i];
			}

			if constexpr (newDim == dimension + 1) {
				result[dimension] = w;
			}

			return result;
		}

		// 从齐次坐标转换回普通坐标
		template<size_t newDim = dimension - 1>
		Vector<T, newDim> fromHomogeneous() const
		{
			static_assert(newDim < dimension, "New dimension must be smaller");
			static_assert(dimension > 1, "Cannot reduce dimension of vector less than 2 dimensions");

			Vector<T, newDim> result;
			T w = data[dimension - 1];

			// 处理w=0的情况（无穷远点或者方向向量）
			if (std::abs(w) < T{ NEARZERO_THRESHOLD }) {
				// 返回方向向量（不进行透视除法）
				for (size_t i = 0; i < newDim; ++i) {
					result[i] = data[i];
				}
			}
			else {
				// 进行透视除法
				for (size_t i = 0; i < newDim; ++i) 
				{
					result[i] = data[i] / w;
				}
			}

			return result;
		}

		template<typename = std::enable_if_t<dimension == 3>>
		Matrix<T, 4, 4, std::enable_if_t<std::is_arithmetic_v<T>>> toTranslationMatrix(bool usedWithOrient = false)
		{
			Matrix<T, 4, 4, std::enable_if_t<std::is_arithmetic_v<T>>> result = MTXIdentity<T, 4>();
			result(0, 3) = data[0];
			result(1, 3) = data[1];
			result(2, 3) = data[2];

			if (usedWithOrient)
			{
				result(3, 3) = 0;
			}

			return result;
		}

		Matrix<T, dimension, dimension, std::enable_if_t<std::is_arithmetic_v<T>>> skewSymmetric() const
		{

			static_assert(dimension > 1, "dimension must greater than 1");

			if constexpr (dimension == 2) 
			{
				Matrix<T, 2, 2, std::enable_if_t<std::is_arithmetic_v<T>>> result{};
				
				result(0, 1) = -T{1};
				result(1, 0) = T{1};
				
				return result;
			}
			
			Matrix<T, dimension, dimension, std::enable_if_t<std::is_arithmetic_v<T>>> result{};  // 初始化为零矩阵

			if constexpr (dimension == 3) 
			{
				// 三维向量的标准反对矩阵
				result(0, 1) = -data[2];  // -z
				result(0, 2) = data[1];   //  y
				result(1, 0) = data[2];   //  z
				result(1, 2) = -data[0];  // -x
				result(2, 0) = -data[1];  // -y
				result(2, 1) = data[0];   //  x
			}
			else 
			{
				// 高维向量的通用构造方法
				// 使用循环填充模式
				size_t idx = 0;
				for (size_t i = 0; i < dimension && i < dimension; ++i) 
				{
					for (size_t j = i + 1; j < dimension && j < dimension; ++j)
					{
						if (idx < dimension) 
						{
							T value = (idx % 2 == 0) ? data[idx] : -data[idx];
							result(i, j) = value;
							result(j, i) = -value;
							idx++;
						}
					}
				}
			}

			return result;
		}

		template<size_t D = dimension>
		typename std::enable_if_t<D == 2, Matrix<T, 3, 3, std::enable_if_t<std::is_arithmetic_v<T>>>> skewSymmetric_2D() const
		{

			static_assert(dimension == 2, "dimension must be 2");
			
			Matrix<T, 3, 3, std::enable_if_t<std::is_arithmetic_v<T>>> result{};  // 初始化为零矩阵
			
			result(0, 1) = T{ 0 };      // 0 (因为z=0)
			result(0, 2) = data[1];   //  y
			result(1, 0) = T{ 0 };      //  0 (因为z=0)
			result(1, 2) = -data[0];  // -x
			result(2, 0) = -data[1];  // -y
			result(2, 1) = data[0];   //  x

			return result;
			
			
		}
		
/*
* 成员函数：复合赋值运算符
*
*/
		virtual Vector<T, dimension>& operator=(const Vector<T, dimension>& other)
		{
			for (size_t idx = 0; idx < dimension; ++idx)
			{
				this->data[idx] = other[idx];
			}

			return *this;
		}

		virtual Vector<T, dimension>& operator=(Vector<T, dimension>&& other)
		{
			for (size_t idx = 0; idx < dimension; ++idx)
			{
				this->data[idx] = other[idx];
			}

			return *this;
		}
		
		virtual Vector<T, dimension>& operator+=(const Vector<T, dimension>& other)
		{
			for (size_t i = 0; i < dimension; ++i) {
				data[i] += other.data[i];
			}
			return *this;
		}

		virtual Vector<T, dimension>& operator+=(const T& scalar)
		{
			for (size_t i = 0; i < dimension; ++i) {
				data[i] += scalar;
			}
			return *this;
		}

		virtual Vector<T, dimension>& operator-=(const Vector<T, dimension>& other)
		{
			for (size_t i = 0; i < dimension; ++i) {
				data[i] -= other.data[i];
			}
			return *this;
		}

		virtual Vector<T, dimension>& operator-=(const T& scalar)
		{
			for (size_t i = 0; i < dimension; ++i) {
				data[i] -= scalar;
			}
			return *this;
		}

		virtual Vector<T, dimension>& operator*=(const Vector<T, dimension>& other)
		{
			for (size_t i = 0; i < dimension; ++i) {
				data[i] *= other.data[i];
			}
			return *this;
		}

		virtual Vector<T, dimension>& operator*=(const T& scalar)
		{
			for (size_t i = 0; i < dimension; ++i) {
				data[i] *= scalar;
			}
			return *this;
		}

		virtual Vector<T, dimension>& operator/=(const Vector<T, dimension>& other)
		{
			// assert 移到循环外，避免循环内分支阻碍编译器向量化
			for (size_t i = 0; i < dimension; ++i)
				assert(other.data[i] != T{ 0 } && "Division by zero");
			for (size_t i = 0; i < dimension; ++i)
				data[i] /= other.data[i];
			return *this;
		}

		virtual Vector<T, dimension>& operator/=(const T& scalar)
		{
			assert(scalar != T{ 0 } && "Division by zero");
			for (size_t i = 0; i < dimension; ++i) {
				data[i] /= scalar;
			}
			return *this;
		}

	private:

/*
* 友元函数：二元运算符
*
*/
		friend Vector<T, dimension> operator+(const Vector<T, dimension>& lhs, const Vector<T, dimension>& rhs)
		{
			Vector result;
			for (size_t i = 0; i < dimension; ++i)
			{
				result.data[i] = lhs.data[i] + rhs.data[i];
			}
			return result;
		}

		friend Vector<T, dimension> operator+(const Vector<T, dimension>& vec, const T& scalar)
		{
			Vector result;
			for (size_t i = 0; i < dimension; ++i) {
				result.data[i] = vec.data[i] + scalar;
			}
			return result;
		}

		friend Vector<T, dimension> operator+(const T& scalar, const Vector<T, dimension>& vec)
		{
			return vec + scalar;
		}

		friend Vector<T, dimension> operator-(const Vector<T, dimension>& lhs, const Vector<T, dimension>& rhs)
		{
			Vector<T, dimension> result;
			for (size_t i = 0; i < dimension; ++i)
			{
				result.data[i] = lhs.data[i] - rhs.data[i];
			}
			return result;
		}

		friend Vector<T, dimension> operator-(const Vector<T, dimension>& vec, const T& scalar)
		{
			Vector<T, dimension> result;
			for (size_t i = 0; i < dimension; ++i)
			{
				result.data[i] = vec.data[i] - scalar;
			}
			return result;
		}

		friend Vector<T, dimension> operator-(const Vector<T, dimension>& vec)
		{
			Vector<T, dimension> result;
			for (size_t i = 0; i < dimension; ++i)
			{
				result.data[i] = -vec.data[i];
			}
			return result;
		}

		friend Vector<T, dimension> operator*(const Vector<T, dimension>& lhs, const Vector<T, dimension>& rhs)
		{
			Vector<T, dimension> result;
			for (size_t i = 0; i < dimension; ++i) 
			{
				result.data[i] = lhs.data[i] * rhs.data[i];
			}
			return result;
		}

		friend Vector<T, dimension> operator*(const Vector<T, dimension>& vec, const T& scalar)
		{
			Vector<T, dimension> result;
			for (size_t i = 0; i < dimension; ++i) 
			{
				result.data[i] = vec.data[i] * scalar;
			}
			return result;
		}

		friend Vector<T, dimension> operator*(const T& scalar, const Vector<T, dimension>& vec)
		{
			return vec * scalar;
		}

		friend Vector<T, dimension> operator/(const Vector<T, dimension>& lhs, const Vector<T, dimension>& rhs)
		{
			Vector<T, dimension> result;
			for (size_t i = 0; i < dimension; ++i) 
			{
				assert(rhs.data[i] != T{ 0 } && "Division by zero");
				result.data[i] = lhs.data[i] / rhs.data[i];
			}
			return result;
		}

		friend Vector<T, dimension> operator/(const Vector<T, dimension>& vec, const T& scalar)
		{
			assert(scalar != T{ 0 } && "Division by zero");
			Vector<T, dimension> result;
			for (size_t i = 0; i < dimension; ++i) 
			{
				result.data[i] = vec.data[i] / scalar;
			}
			return result;
		}

		// 比较运算符
		friend bool operator==(const Vector<T, dimension>& lhs, const Vector<T, dimension>& rhs)
		{
			if constexpr (std::is_floating_point_v<T>) 
			{
				// 浮点数使用epsilon比较
				for (size_t i = 0; i < dimension; ++i)
				{
					if (std::abs(lhs.data[i] - rhs.data[i]) > T{ NEARZERO_THRESHOLD })
					{
						return false;
					}
				}
			}
			else 
			{
				// 整数类型直接比较
				for (size_t i = 0; i < dimension; ++i)
				{
					if (lhs.data[i] != rhs.data[i])
					{
						return false;
					}
				}
			}

			return true;
		}

		friend bool operator!=(const Vector<T, dimension>& lhs, const Vector<T, dimension>& rhs)
		{
			return !(lhs == rhs);
		}

	public:
		union
		{
			alignas(16) std::array<T, dimension> data;
			Swizzle<T, dimension, 0> x;
			Swizzle<T, dimension, 1> y;
			Swizzle<T, dimension, 2> z;
			Swizzle<T, dimension, 3> w;
			Swizzle<T, dimension, 0, 1> xy;
			Swizzle<T, dimension, 1, 0> yx;
			Swizzle<T, dimension, 0, 2> xz;
			Swizzle<T, dimension, 2, 0> zx;
			Swizzle<T, dimension, 0, 3> xw;
			Swizzle<T, dimension, 3, 0> wx;
			Swizzle<T, dimension, 1, 2> yz;
			Swizzle<T, dimension, 2, 1> zy;
			Swizzle<T, dimension, 1, 3> yw;
			Swizzle<T, dimension, 3, 1> wy;
			Swizzle<T, dimension, 2, 3> zw;
			Swizzle<T, dimension, 3, 2> wz;
			Swizzle<T, dimension, 0, 1, 2> xyz;
			Swizzle<T, dimension, 0, 2, 1> xzy;
			Swizzle<T, dimension, 1, 0, 2> yxz;
			Swizzle<T, dimension, 1, 2, 0> yzx;
			Swizzle<T, dimension, 2, 0, 1> zxy;
			Swizzle<T, dimension, 2, 1, 0> zyx;
			Swizzle<T, dimension, 0, 1, 3> xyw;
			Swizzle<T, dimension, 0, 3, 1> xwy;
			Swizzle<T, dimension, 1, 0, 3> yxw;
			Swizzle<T, dimension, 1, 3, 0> ywx;
			Swizzle<T, dimension, 3, 0, 1> wxy;
			Swizzle<T, dimension, 3, 1, 0> wyx;
			Swizzle<T, dimension, 1, 2, 3> yzw;
			Swizzle<T, dimension, 1, 3, 2> ywz;
			Swizzle<T, dimension, 2, 1, 3> zyw;
			Swizzle<T, dimension, 2, 3, 1> zwy;
			Swizzle<T, dimension, 3, 1, 2> wyz;
			Swizzle<T, dimension, 3, 2, 1> wzy;
			Swizzle<T, dimension, 0, 2, 3> xzw;
			Swizzle<T, dimension, 0, 3, 2> xwz;
			Swizzle<T, dimension, 2, 0, 3> zxw;
			Swizzle<T, dimension, 2, 3, 0> zwx;
			Swizzle<T, dimension, 3, 0, 2> wxz;
			Swizzle<T, dimension, 3, 2, 0> wzx;
			Swizzle<T, dimension, 0, 1, 2, 3> xyzw;
			Swizzle<T, dimension, 0, 1, 3, 2> xywz;
			Swizzle<T, dimension, 0, 2, 1, 3> xzyw;
			Swizzle<T, dimension, 0, 2, 3, 1> xzwy;
			Swizzle<T, dimension, 0, 3, 1, 2> xwyz;
			Swizzle<T, dimension, 0, 3, 2, 1> xwzy;
			Swizzle<T, dimension, 1, 0, 2, 3> yxzw;
			Swizzle<T, dimension, 1, 0, 3, 2> yxwz;
			Swizzle<T, dimension, 1, 2, 0, 3> yzxw;
			Swizzle<T, dimension, 1, 2, 3, 0> yzwx;
			Swizzle<T, dimension, 1, 3, 0, 2> ywxz;
			Swizzle<T, dimension, 1, 3, 2, 0> ywzx;
			Swizzle<T, dimension, 2, 0, 1, 3> zxyw;
			Swizzle<T, dimension, 2, 0, 3, 1> zxwy;
			Swizzle<T, dimension, 2, 1, 0, 3> zyxw;
			Swizzle<T, dimension, 2, 1, 3, 0> zywx;
			Swizzle<T, dimension, 2, 3, 0, 1> zwxy;
			Swizzle<T, dimension, 2, 3, 1, 0> zwyx;
			Swizzle<T, dimension, 3, 0, 1, 2> wxyz;
			Swizzle<T, dimension, 3, 0, 2, 1> wxzy;
			Swizzle<T, dimension, 3, 1, 0, 2> wyxz;
			Swizzle<T, dimension, 3, 1, 2, 0> wyzx;
			Swizzle<T, dimension, 3, 2, 0, 1> wzxy;
			Swizzle<T, dimension, 3, 2, 1, 0> wzyx;

			Swizzle<T, dimension, 0> r;
			Swizzle<T, dimension, 1> g;
			Swizzle<T, dimension, 2> b;
			Swizzle<T, dimension, 3> a;
			Swizzle<T, dimension, 0, 1> rg;
			Swizzle<T, dimension, 1, 0> gr;
			Swizzle<T, dimension, 0, 2> rb;
			Swizzle<T, dimension, 2, 0> br;
			Swizzle<T, dimension, 0, 3> ra;
			Swizzle<T, dimension, 3, 0> ar;
			Swizzle<T, dimension, 1, 2> gb;
			Swizzle<T, dimension, 2, 1> bg;
			Swizzle<T, dimension, 1, 3> ga;
			Swizzle<T, dimension, 3, 1> ag;
			Swizzle<T, dimension, 2, 3> ba;
			Swizzle<T, dimension, 3, 2> ab;
			Swizzle<T, dimension, 0, 1, 2> rgb;
			Swizzle<T, dimension, 0, 2, 1> rbg;
			Swizzle<T, dimension, 1, 0, 2> grb;
			Swizzle<T, dimension, 1, 2, 0> gbr;
			Swizzle<T, dimension, 2, 0, 1> brg;
			Swizzle<T, dimension, 2, 1, 0> bgr;
			Swizzle<T, dimension, 0, 1, 3> rga;
			Swizzle<T, dimension, 0, 3, 1> rag;
			Swizzle<T, dimension, 1, 0, 3> gra;
			Swizzle<T, dimension, 1, 3, 0> gar;
			Swizzle<T, dimension, 3, 0, 1> arg;
			Swizzle<T, dimension, 3, 1, 0> agr;
			Swizzle<T, dimension, 1, 2, 3> gba;
			Swizzle<T, dimension, 1, 3, 2> gab;
			Swizzle<T, dimension, 2, 1, 3> bga;
			Swizzle<T, dimension, 2, 3, 1> bag;
			Swizzle<T, dimension, 3, 1, 2> agb;
			Swizzle<T, dimension, 3, 2, 1> abg;
			Swizzle<T, dimension, 0, 2, 3> rba;
			Swizzle<T, dimension, 0, 3, 2> rab;
			Swizzle<T, dimension, 2, 0, 3> bra;
			Swizzle<T, dimension, 2, 3, 0> bar;
			Swizzle<T, dimension, 3, 0, 2> arb;
			Swizzle<T, dimension, 3, 2, 0> abr;
			Swizzle<T, dimension, 0, 1, 2, 3> rgba;
			Swizzle<T, dimension, 0, 1, 3, 2> rgab;
			Swizzle<T, dimension, 0, 2, 1, 3> rbga;
			Swizzle<T, dimension, 0, 2, 3, 1> rbag;
			Swizzle<T, dimension, 0, 3, 1, 2> ragb;
			Swizzle<T, dimension, 0, 3, 2, 1> rabg;
			Swizzle<T, dimension, 1, 0, 2, 3> grba;
			Swizzle<T, dimension, 1, 0, 3, 2> grab;
			Swizzle<T, dimension, 1, 2, 0, 3> gbra;
			Swizzle<T, dimension, 1, 2, 3, 0> gbar;
			Swizzle<T, dimension, 1, 3, 0, 2> garb;
			Swizzle<T, dimension, 1, 3, 2, 0> gabr;
			Swizzle<T, dimension, 2, 0, 1, 3> brga;
			Swizzle<T, dimension, 2, 0, 3, 1> brag;
			Swizzle<T, dimension, 2, 1, 0, 3> bgra;
			Swizzle<T, dimension, 2, 1, 3, 0> bgar;
			Swizzle<T, dimension, 2, 3, 0, 1> barg;
			Swizzle<T, dimension, 2, 3, 1, 0> bagr;
			Swizzle<T, dimension, 3, 0, 1, 2> argb;
			Swizzle<T, dimension, 3, 0, 2, 1> arbg;
			Swizzle<T, dimension, 3, 1, 0, 2> agrb;
			Swizzle<T, dimension, 3, 1, 2, 0> agbr;
			Swizzle<T, dimension, 3, 2, 0, 1> abrg;
			Swizzle<T, dimension, 3, 2, 1, 0> abgr;
		};
		
	};

/*
* 按照标准化的方式在控制台输出Vector
* 其直接调用Vector中内置的print函数
*/
	template<typename T, size_t dimension>
	std::ostream& operator<<(std::ostream& out, const Vector<T,dimension>& vec)
	{
		return vec.print(out);
	}

/*
 * 全局函数：向量运算
 * 
 */

	template<typename T, size_t D>
	T dot(const Vector<T, D>& a, const Vector<T, D>& b)
	{
		T result = T{ 0 };
		for (size_t i = 0; i < D; ++i) 
		{
			result += a[i] * b[i];
		}
		return result;
	}

	template<typename T, size_t Dim1, size_t Dim2, int... Indices1, int... Indices2, size_t dim_a = sizeof...(Indices1), size_t dim_b = sizeof...(Indices2)>
	std::enable_if_t<dim_a == dim_b,T> dot(const Swizzle<T, Dim1, Indices1...>& a, const Swizzle<T, Dim2, Indices2...>& b)
	{
		if constexpr (dim_a > 1)
		{
			Vector<T, dim_a> va = a;
			Vector<T, dim_b> vb = b;
			return dot(va, vb);
		}
		else
		{
			T va = a;
			T vb = b;
			return va * vb;
		}
	}

	// 3D叉积
	template<typename T, size_t dimension>
	typename std::enable_if_t<dimension == 3, Vector<T, dimension>>
	cross(const Vector<T, dimension>& a, const Vector<T, dimension>& b)
	{
		return {
			a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x
		};
	}

	template<typename T, size_t Dim1, size_t Dim2, int... Indices1, int... Indices2, size_t dim_a = sizeof...(Indices1), size_t dim_b = sizeof...(Indices2)>
	std::enable_if_t<dim_a == 3 && dim_b == 3, Vector<T, 3>> cross(const Swizzle<T, Dim1, Indices1...>& a, const Swizzle<T, Dim2, Indices2...>& b)
	{
	
		Vector<T, 3> va = a;
		Vector<T, 3> vb = b;
		return cross(va, vb);
		
	}

	// 2D叉积（返回标量）
	template<typename T, size_t dimension>
	typename std::enable_if_t<dimension == 2, T>
	cross(const Vector<T, dimension>& a, const Vector<T, dimension>& b)
	{
		return a.x * b.y - a.y * b.x;
	}

	template<typename T, size_t Dim1, size_t Dim2, int... Indices1, int... Indices2, size_t dim_a = sizeof...(Indices1), size_t dim_b = sizeof...(Indices2)>
	std::enable_if_t<dim_a == 2 && dim_b == 2, T> cross(const Swizzle<T, Dim1, Indices1...>& a, const Swizzle<T, Dim2, Indices2...>& b)
	{
	
		Vector<T, 2> va = a;
		Vector<T, 2> vb = b;
		return cross(va, vb);
		
	}

	// 距离函数
	template<typename T, size_t D>
	T distance(const Vector<T, D>& a, const Vector<T, D>& b, bool dimensionalityReduction  = true)
	{
		return (b - a).length(dimensionalityReduction);
	}

	template<typename T, size_t Dim1, size_t Dim2, int... Indices1, int... Indices2, size_t dim_a = sizeof...(Indices1), size_t dim_b = sizeof...(Indices2)>
	std::enable_if_t<dim_a == dim_b, T> distance(const Swizzle<T, Dim1, Indices1...>& a, const Swizzle<T, Dim2, Indices2...>& b, bool dimensionalityReduction  = true)
	{
		Vector<T, dim_a> va = a;
		Vector<T, dim_b> vb = b;
		return distance(va, vb, dimensionalityReduction);
	}

	template<typename T, size_t D>
	T distanceSquared(const Vector<T, D>& a, const Vector<T, D>& b, bool dimensionalityReduction = true)
	{
		return (b - a).lengthSquared(dimensionalityReduction);
	}

	template<typename T, size_t Dim1, size_t Dim2, int... Indices1, int... Indices2, size_t dim_a = sizeof...(Indices1), size_t dim_b = sizeof...(Indices2)>
	std::enable_if_t<dim_a == dim_b, T> distanceSquared(const Swizzle<T, Dim1, Indices1...>& a, const Swizzle<T, Dim2, Indices2...>& b, bool dimensionalityReduction  = true)
	{
		Vector<T, dim_a> va = a;
		Vector<T, dim_b> vb = b;
		return distanceSquared(va, vb, dimensionalityReduction);
	}

	// 角度函数
	template<typename T, size_t D>
	T angle(const Vector<T, D>& a, const Vector<T, D>& b, bool dimensionalityReduction = true)
	{
		T dot_product = dot(a, b);
		T magnitude_product = a.length(dimensionalityReduction) * b.length(dimensionalityReduction);
		if (magnitude_product == T{ 0 }) {
			return T{ 0 };
		}
		return std::acos(std::clamp(dot_product / magnitude_product, T{ -1 }, T{ 1 }));
	}

	template<typename T, size_t Dim1, size_t Dim2, int... Indices1, int... Indices2, size_t dim_a = sizeof...(Indices1), size_t dim_b = sizeof...(Indices2)>
	std::enable_if_t<dim_a == dim_b, T> angle(const Swizzle<T, Dim1, Indices1...>& a, const Swizzle<T, Dim2, Indices2...>& b, bool dimensionalityReduction  = true)
	{
		Vector<T, dim_a> va = a;
		Vector<T, dim_b> vb = b;
		return angle(va,vb,dimensionalityReduction);
	}

	// 线性插值
	template<typename T, size_t D>
	Vector<T, D> lerp(const Vector<T, D>& a, const Vector<T, D>& b, T t)
	{
		return a.lerp(b, t);
	}

	template<typename T, size_t Dim1, size_t Dim2, int... Indices1, int... Indices2, size_t dim_a = sizeof...(Indices1), size_t dim_b = sizeof...(Indices2)>
	std::enable_if_t<dim_a == dim_b, T> lerp(const Swizzle<T, Dim1, Indices1...>& a, const Swizzle<T, Dim2, Indices2...>& b, T t)
	{
		Vector<T, dim_a> va = a;
		Vector<T, dim_b> vb = b;
		return lerp(va,vb,t);
	}

	// 球面线性插值（仅用于单位向量）
	template<typename T, size_t D>
	Vector<T, D> slerp(const Vector<T, D>& a, const Vector<T, D>& b, T t)
	{
		T dot_product = std::clamp(dot(a, b), T{ -1 }, T{ 1 });
		T theta = std::acos(std::abs(dot_product));

		if (theta < T{ 1e-6 }) {
			return lerp(a, b, t);  // 向量几乎平行，使用线性插值
		}

		T sin_theta = std::sin(theta);
		T factor_a = std::sin((T{ 1 } - t) * theta) / sin_theta;
		T factor_b = std::sin(t * theta) / sin_theta;

		if (dot_product < T{ 0 }) {
			factor_b = -factor_b;  // 选择较短路径
		}

		return factor_a * a + factor_b * b;
	}

	template<typename T, size_t Dim1, size_t Dim2, int... Indices1, int... Indices2, size_t dim_a = sizeof...(Indices1), size_t dim_b = sizeof...(Indices2)>
	std::enable_if_t<dim_a == dim_b, T> slerp(const Swizzle<T, Dim1, Indices1...>& a, const Swizzle<T, Dim2, Indices2...>& b, T t)
	{
		Vector<T, dim_a> va = a;
		Vector<T, dim_b> vb = b;
		return slerp(va,vb,t);
	}

	
/*
 * Swizzle
 * 
 */
	
	template <typename T, size_t dimension, int... Indices>
	struct Swizzle
	{
	    // 编译期计算 swizzle 的维度
	    static constexpr size_t swizzle_dim = sizeof...(Indices);
	    
	    // 编译期检查所有索引是否在范围内
	    static constexpr bool check_indices()
	    {
	    	return ((Indices < static_cast<int>(dimension)) && ...);
	    }
	    
	public:
	    // ==================== 1D Swizzle 特化功能 ====================
	    
	    // 1D: 支持赋值单个标量
	    template<size_t D = swizzle_dim>
	    std::enable_if_t<D == 1, Swizzle&> operator=(const T& value)
	    {
	        static_assert(check_indices(), 
	            "Swizzle index out of bounds! "
	            "You're trying to access a component that doesn't exist in this vector.");
	        
	        constexpr int idx_array[] = {Indices...};
	        elem(idx_array[0]) = value;
	        return *this;
	    }
	    
	    // 1D: 隐式转换为 T
	    template<size_t D = swizzle_dim>
	    operator std::enable_if_t<D == 1, T>() const
	    {
	        static_assert(check_indices(), 
	            "Swizzle index out of bounds! "
	            "You're trying to access a component that doesn't exist in this vector.");
	        
	        constexpr int idx_array[] = {Indices...};
	        return elem(idx_array[0]);
	    }
	    
	    // ==================== 通用功能（所有维度）====================
	    
	    // 赋值操作符：从 Vector 赋值
		template<size_t D = swizzle_dim>
		std::enable_if_t<D >= 2, Swizzle&> operator=(const Vector<T, swizzle_dim>& vec)
	    {
	        static_assert(check_indices(),
	            "Swizzle index out of bounds! "
	            "You're trying to access a component that doesn't exist in this vector.");
	        
	        assign_vector_impl(vec, std::make_index_sequence<swizzle_dim>{});
	        return *this;
	    }
	    
	    // 转换为 Vector
		template<size_t D = swizzle_dim>
		operator std::enable_if_t<D >= 2, Vector<T, swizzle_dim>> () const
	    {
	        static_assert(check_indices(), 
	            "Swizzle index out of bounds! "
	            "You're trying to access a component that doesn't exist in this vector.");
	        
	        return make_vector_impl(std::make_index_sequence<swizzle_dim>{});
	    }
	    
	    // 输出操作符
	    friend std::ostream& operator<<(std::ostream& out, const Swizzle& sw)
	    {
	        static_assert(check_indices(), 
	            "Swizzle index out of bounds! "
	            "You're trying to access a component that doesn't exist in this vector.");
	        
	        if constexpr (swizzle_dim == 1)
	        {
	            return out << T(sw);
	        }
	        else
	        {
	            return out << Vector<T, swizzle_dim>(sw);
	        }
	    }

		// ==================== 算术运算符重载 ====================

		// ---------------- 1D Swizzle 与 1D Swizzle 运算（返回标量） ----------------
		
		template<size_t D = swizzle_dim, size_t otherDimension, int OtherIndex>
		typename std::enable_if_t<D == 1, T> 
		operator+(const Swizzle<T, otherDimension, OtherIndex>& other) const
	    {
	    	static_assert(check_indices(), "left swizzle index out of bounds");
	    	static_assert(other.check_indices(), "right swizzle index out of bounds");
	    	return T(*this) + T(other);
	    }
		
		template<size_t D = swizzle_dim, size_t otherDimension, int OtherIndex>
		typename std::enable_if_t<D == 1, T> 
		operator-(const Swizzle<T, otherDimension, OtherIndex>& other) const
	    {
	    	static_assert(check_indices(), "left swizzle index out of bounds");
	    	static_assert(other.check_indices(), "right swizzle index out of bounds");
	    	return T(*this) - T(other);
	    }
		
		template<size_t D = swizzle_dim, size_t otherDimension, int OtherIndex>
		typename std::enable_if_t<D == 1, T> 
		operator*(const Swizzle<T, otherDimension, OtherIndex>& other) const
	    {
	    	static_assert(check_indices(), "left swizzle index out of bounds");
	    	static_assert(other.check_indices(), "right swizzle index out of bounds");
	    	return T(*this) * T(other);
	    }
		
		template<size_t D = swizzle_dim, size_t otherDimension, int OtherIndex>
		typename std::enable_if_t<D == 1, T> 
		operator/(const Swizzle<T, otherDimension, OtherIndex>& other) const
	    {
	    	static_assert(check_indices(), "left swizzle index out of bounds");
	    	static_assert(other.check_indices(), "right swizzle index out of bounds");
	    	return T(*this) / T(other);
	    }
		
		// ---------------- 1D Swizzle 与不同维度 Swizzle 运算（广播） ----------------
		template<size_t D = swizzle_dim, int... OtherIndices>
		typename std::enable_if_t<D == 1 && (sizeof...(OtherIndices) > 1), Vector<T, sizeof...(OtherIndices)>> 
		operator+(const Swizzle<T, dimension, OtherIndices...>& other) const
	    {
	    	static_assert(check_indices(), "left swizzle index out of bounds");
	    	static_assert(other.check_indices(), "right swizzle index out of bounds");
	    	constexpr size_t other_dim = sizeof...(OtherIndices);
	    	return Vector<T, other_dim>(other) + T(*this);
	    }
		
		template<size_t D = swizzle_dim, int... OtherIndices>
		typename std::enable_if_t<D == 1 && (sizeof...(OtherIndices) > 1), Vector<T, sizeof...(OtherIndices)>> 
		operator*(const Swizzle<T, dimension, OtherIndices...>& other) const
	    {
	    	static_assert(check_indices(), "left swizzle index out of bounds");
	    	static_assert(other.check_indices(), "right swizzle index out of bounds");
	    	constexpr size_t other_dim = sizeof...(OtherIndices);
	    	return Vector<T, other_dim>(other) * T(*this);
	    }
		
		// ---------------- Swizzle 与 Swizzle 运算（相同IndicesNumber） ----------------

		// swizzle + 1D swizzle（右侧1D广播）

		template<int OtherIndex, size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> 
		operator+(const Swizzle<T, dimension, OtherIndex>& other) const
	    {
	    	static_assert(check_indices(), "Left swizzle index out of bounds");
	        
	    	Vector<T, swizzle_dim> lhs = *this;
	    	T scalar_value = T(other);
	    	return lhs + scalar_value;
	    }

		template<int OtherIndex, size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> 
		operator-(const Swizzle<T, dimension, OtherIndex>& other) const
	    {
	    	static_assert(check_indices(), "Left swizzle index out of bounds");
	        
	    	Vector<T, swizzle_dim> lhs = *this;
	    	T scalar_value = T(other);
	    	return lhs - scalar_value;
	    }

		template<int OtherIndex, size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> 
		operator*(const Swizzle<T, dimension, OtherIndex>& other) const
	    {
	    	static_assert(check_indices(), "Left swizzle index out of bounds");
	        
	    	Vector<T, swizzle_dim> lhs = *this;
	    	T scalar_value = T(other);
	    	return lhs * scalar_value;
	    }

		template<int OtherIndex, size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> 
		operator/(const Swizzle<T, dimension, OtherIndex>& other) const
	    {
	    	static_assert(check_indices(), "Left swizzle index out of bounds");
	        
	    	Vector<T, swizzle_dim> lhs = *this;
	    	T scalar_value = T(other);
	    	return lhs / scalar_value;
	    }

		// Swizzle + Swizzle

		template<int ... OtherIndices, size_t OtherDimension, size_t D = swizzle_dim, size_t OtherD = sizeof ... (OtherIndices)>
		typename std::enable_if_t<(D > 1 && D == OtherD), Vector<T, swizzle_dim>> 
		operator+(const Swizzle<T, OtherDimension, OtherIndices ...>& other) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	Vector<T, swizzle_dim> lhs = *this;
	    	Vector<T, swizzle_dim> rhs = other;
	    	return lhs + rhs;
	    }

		template<int ... OtherIndices, size_t OtherDimension, size_t D = swizzle_dim, size_t OtherD = sizeof ... (OtherIndices)>
		typename std::enable_if_t<(D > 1 && D == OtherD), Vector<T, swizzle_dim>> 
		operator-(const Swizzle<T, OtherDimension, OtherIndices ... >& other) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	Vector<T, swizzle_dim> lhs = *this;
	    	Vector<T, swizzle_dim> rhs = other;
	    	return lhs - rhs;
	    }

		template<int ... OtherIndices, size_t OtherDimension, size_t D = swizzle_dim, size_t OtherD = sizeof ... (OtherIndices)>
		typename std::enable_if_t<(D > 1 && D == OtherD), Vector<T, swizzle_dim>> 
		operator*(const Swizzle<T, OtherDimension, OtherIndices ...>& other) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	Vector<T, swizzle_dim> lhs = *this;
	    	Vector<T, swizzle_dim> rhs = other;
	    	return lhs * rhs;
	    }

		template<int ... OtherIndices, size_t OtherDimension, size_t D = swizzle_dim, size_t OtherD = sizeof ... (OtherIndices)>
		typename std::enable_if_t<(D > 1 && D == OtherD), Vector<T, swizzle_dim>> 
		operator/(const Swizzle<T, OtherDimension, OtherIndices ... >& other) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	Vector<T, swizzle_dim> lhs = *this;
	    	Vector<T, swizzle_dim> rhs = other;
	    	return lhs / rhs;
	    }

		// ---------------- 1D Swizzle 与 Vector 运算 ----------------
		
		template<size_t D = swizzle_dim, size_t VDimension>
		typename std::enable_if_t<D == 1, Vector<T,VDimension>> operator+(const Vector<T,VDimension>& vec) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec + T(*this);
	    }

		template<size_t D = swizzle_dim, size_t VDimension>
		typename std::enable_if_t<D == 1, Vector<T,VDimension>> operator*(const Vector<T,VDimension>& vec) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec * T(*this);
	    }

		// ---------------- Swizzle 与 Vector 运算 ----------------

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T,swizzle_dim>> operator+(const Vector<T,swizzle_dim>& vec) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec + Vector<T,swizzle_dim>(*this);
	    }

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T,swizzle_dim>> operator-(const Vector<T,swizzle_dim>& vec) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return Vector<T,swizzle_dim>(*this) - vec;
	    }

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T,swizzle_dim>> operator*(const Vector<T,swizzle_dim>& vec) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec * Vector<T,swizzle_dim>(*this);
	    }

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T,swizzle_dim>> operator/(const Vector<T,swizzle_dim>& vec) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return Vector<T,swizzle_dim>(*this) / vec;
	    }

		// ---------------- 1D Swizzle 与 Scalar 运算 ----------------

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<D == 1, T> operator+(const T& scalar) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return T(*this) + scalar;
	    }

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<D == 1, T> operator-(const T& scalar) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return T(*this) - scalar;
	    }

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<D == 1, T> operator*(const T& scalar) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return T(*this) * scalar;
	    }

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<D == 1, T> operator/(const T& scalar) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return T(*this) / scalar;
	    }
		
		// ---------------- Swizzle 与 Scalar 运算 ----------------
		
		template<size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> operator+(const T& scalar) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	Vector<T, swizzle_dim> lhs = *this;
	    	return lhs + scalar;
	    }
		
		template<size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> operator-(const T& scalar) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	Vector<T, swizzle_dim> lhs = *this;
	    	return lhs - scalar;
	    }

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> operator*(const T& scalar) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	Vector<T, swizzle_dim> lhs = *this;
	    	return lhs * scalar;
	    }

		template<size_t D = swizzle_dim>
		typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> operator/(const T& scalar) const
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	Vector<T, swizzle_dim> lhs = *this;
	    	return lhs / scalar;
	    }

		// ---------------- 友元函数：支持反向运算 ----------------
		
		// Vector & Swizzle
		template<size_t D = swizzle_dim>
		friend typename std::enable_if_t<(D > 1), Vector<T,swizzle_dim>> operator+(const Vector<T,swizzle_dim>& vec, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec + Vector<T,swizzle_dim>(sw);
	    }

		template<size_t D = swizzle_dim>
		friend typename std::enable_if_t<(D > 1), Vector<T,swizzle_dim>> operator-(const Vector<T,swizzle_dim>& vec, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec - Vector<T,swizzle_dim>(sw);
	    }

		template<size_t D = swizzle_dim>
		friend typename std::enable_if_t<(D > 1), Vector<T,swizzle_dim>> operator*(const Vector<T,swizzle_dim>& vec, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec * Vector<T,swizzle_dim>(sw);
	    }

		template<size_t D = swizzle_dim>
		friend typename std::enable_if_t<(D > 1), Vector<T,swizzle_dim>> operator/(const Vector<T,swizzle_dim>& vec, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec / Vector<T,swizzle_dim>(sw) ;
	    }
		
		// Vector & 1D Swizzle（广播）
		template<size_t D = swizzle_dim, size_t VDimension>
		friend typename std::enable_if_t<(D > 1), Vector<T,VDimension>> operator+(const Vector<T,VDimension>& vec, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return sw + vec;
	    }

		template<size_t D = swizzle_dim, size_t VDimension>
		friend typename std::enable_if_t<(D > 1), Vector<T,VDimension>> operator-(const Vector<T,VDimension>& vec, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec - T(sw);
	    }

		template<size_t D = swizzle_dim, size_t VDimension>
		friend typename std::enable_if_t<(D > 1), Vector<T,VDimension>> operator*(const Vector<T,VDimension>& vec, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return sw * vec;
	    }

		template<size_t D = swizzle_dim, size_t VDimension>
		friend typename std::enable_if_t<(D > 1), Vector<T,VDimension>> operator/(const Vector<T,VDimension>& vec, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return vec / T(sw);
	    }
		
		// Scalar & Swizzle

		template<size_t D = swizzle_dim>
		friend  typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>> operator+(const T& scalar, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return sw + scalar;
	    }
		
		template<size_t D = swizzle_dim>
		friend  typename std::enable_if_t<(D > 1), Vector<T, swizzle_dim>>  operator*(const T& scalar, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return sw * scalar;
	    }
		
		// Scalar & 1D Swizzle
		
		template<size_t D = swizzle_dim>
		friend typename std::enable_if_t<D == 1, T> operator+(const T& scalar, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return sw + scalar;
	    }
		
		template<size_t D = swizzle_dim>
		friend typename std::enable_if_t<D == 1, T> operator-(const T& scalar, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return scalar - T(sw) ;
	    }

		template<size_t D = swizzle_dim>
		friend typename std::enable_if_t<D == 1, T>  operator*(const T& scalar, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return sw * scalar;
	    }

		template<size_t D = swizzle_dim>
		friend typename std::enable_if_t<D == 1, T> operator/(const T& scalar, const Swizzle& sw)
	    {
	    	static_assert(check_indices(), "Swizzle index out of bounds");
	    	return scalar/T(sw);
	    }
		

	protected:
		// 访问底层数据的辅助函数
		// 注意：这些函数不能是 static 的，因为需要通过 this 指针访问 union 中的数据
		T elem(int i) const
		{
			// 通过 reinterpret_cast 将 this 指针转换为 T* 来访问 union 中的 data 数组
			return reinterpret_cast<const T*>(this)[i];
		}
	    
		T& elem(int i)
		{
			return reinterpret_cast<T*>(this)[i];
		}
	    
		// 辅助：使用索引序列创建 Vector
		template<size_t... Is>
		Vector<T, swizzle_dim> make_vector_impl(std::index_sequence<Is...>) const
		{
			// 将 Indices... 展开为数组，然后用索引序列访问
			constexpr int idx_array[] = {Indices...};
			return Vector<T, swizzle_dim>{elem(idx_array[Is])...};
		}
	    
		// 辅助：使用索引序列赋值
		template<size_t... Is>
		void assign_vector_impl(const Vector<T, swizzle_dim>& vec, std::index_sequence<Is...>)
		{
			constexpr int idx_array[] = {Indices...};
			((elem(idx_array[Is]) = vec[Is]), ...);
		}
	};



/*
 * 常用的向量类型定义
 * 
 */

	template <typename T>
	using Vector2 = Vector<T, 2>;

	template <typename T>
	using Vector3 = Vector<T, 3>;

	template <typename T>
	using Vector4 = Vector<T, 4>;

}
