#pragma once

#include "EasyMathAPI.h"
#include "Common.h"

#include <cmath>
#include <iostream>
#include <array>
#include <cassert>
#include <algorithm>
#include <vector>

namespace EM
{
	template <typename T, size_t rows, size_t cols, typename>
	class Matrix;

	template<typename T, typename>
	class Quaternion;

	template <typename T, size_t dimension, int idx>
	struct Swizzle1D;

	template <typename T, size_t dimension, int E0, int E1>
	struct Swizzle2D;
	
	template <typename T, size_t dimension, int E0, int E1, int E2>
	struct Swizzle3D;

	template <typename T, size_t dimension, int E0, int E1, int E2, int E3>
	struct Swizzle4D;

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
		
		static constexpr size_t size() noexcept { return dimension; }
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
		virtual T length(bool dimensionalityReduction = false) const noexcept
		{
			if (!dimensionalityReduction) {
				T sum = T(0);
				for (size_t idx = 0; idx < dimension; ++idx)
				{
					sum += data[idx] * data[idx];
				}
				return std::sqrt(sum);
			}

			T sum = T(0);
			size_t maxDim = std::min(dimension, size_t(3));  // 只计算前3个维度
			for (size_t idx = 0; idx < maxDim; ++idx)
			{
				sum += data[idx] * data[idx];
			}
			return std::sqrt(sum);
		}

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
		
		virtual void normalize(bool dimensionalityReduction = false)
		{
			*this = normalized(dimensionalityReduction);
		}

		virtual bool isNormalized(T epsilon = T{ NEARZERO_THRESHOLD }, bool dimensionalityReduction = true) const noexcept
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
			static_assert(dimension >= 2, "Reflection requires at least 2D vector");

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
			static_assert(dimension > 0, "Cannot reduce dimension of empty vector");

			Vector<T, newDim> result;
			T w = data[dimension - 1];

			// 处理w=0的情况（无穷远点）
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

		virtual Matrix<T, dimension, dimension, std::enable_if_t<std::is_arithmetic_v<T>>> skewSymmetric() const
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
			for (size_t i = 0; i < dimension; ++i) {
				assert(other.data[i] != T{ 0 } && "Division by zero");
				data[i] /= other.data[i];
			}
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
				return true;
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
				return true;
			}
		}

		friend bool operator!=(const Vector<T, dimension>& lhs, const Vector<T, dimension>& rhs)
		{
			return !(lhs == rhs);
		}

	public:
		union
		{
			std::array<T, dimension> data;
			Swizzle1D<T, dimension, 0> x;
			Swizzle1D<T, dimension, 1> y;
			Swizzle1D<T, dimension, 2> z;
			Swizzle1D<T, dimension, 3> w;
			Swizzle2D<T, dimension, 0, 1> xy;
			Swizzle2D<T, dimension, 1, 0> yx;
			Swizzle2D<T, dimension, 0, 2> xz;
			Swizzle2D<T, dimension, 2, 0> zx;
			Swizzle2D<T, dimension, 0, 3> xw;
			Swizzle2D<T, dimension, 3, 0> wx;
			Swizzle2D<T, dimension, 1, 2> yz;
			Swizzle2D<T, dimension, 2, 1> zy;
			Swizzle2D<T, dimension, 1, 3> yw;
			Swizzle2D<T, dimension, 3, 1> wy;
			Swizzle2D<T, dimension, 2, 3> zw;
			Swizzle2D<T, dimension, 3, 2> wz;
			Swizzle3D<T, dimension, 0, 1, 2> xyz;
			Swizzle3D<T, dimension, 0, 2, 1> xzy;
			Swizzle3D<T, dimension, 1, 0, 2> yxz;
			Swizzle3D<T, dimension, 1, 2, 0> yzx;
			Swizzle3D<T, dimension, 2, 0, 1> zxy;
			Swizzle3D<T, dimension, 2, 1, 0> zyx;
			Swizzle3D<T, dimension, 0, 1, 3> xyw;
			Swizzle3D<T, dimension, 0, 3, 1> xwy;
			Swizzle3D<T, dimension, 1, 0, 3> yxw;
			Swizzle3D<T, dimension, 1, 3, 0> ywx;
			Swizzle3D<T, dimension, 3, 0, 1> wxy;
			Swizzle3D<T, dimension, 3, 1, 0> wyx;
			Swizzle3D<T, dimension, 1, 2, 3> yzw;
			Swizzle3D<T, dimension, 1, 3, 2> ywz;
			Swizzle3D<T, dimension, 2, 1, 3> zyw;
			Swizzle3D<T, dimension, 2, 3, 1> zwy;
			Swizzle3D<T, dimension, 3, 1, 2> wyz;
			Swizzle3D<T, dimension, 3, 2, 1> wzy;
			Swizzle3D<T, dimension, 0, 2, 3> xzw;
			Swizzle3D<T, dimension, 0, 3, 2> xwz;
			Swizzle3D<T, dimension, 2, 0, 3> zxw;
			Swizzle3D<T, dimension, 2, 3, 0> zwx;
			Swizzle3D<T, dimension, 3, 0, 2> wxz;
			Swizzle3D<T, dimension, 3, 2, 0> wzx;
			Swizzle4D<T, dimension, 0, 1, 2, 3> xyzw;
			Swizzle4D<T, dimension, 0, 1, 3, 2> xywz;
			Swizzle4D<T, dimension, 0, 2, 1, 3> xzyw;
			Swizzle4D<T, dimension, 0, 2, 3, 1> xzwy;
			Swizzle4D<T, dimension, 0, 3, 1, 2> xwyz;
			Swizzle4D<T, dimension, 0, 3, 2, 1> xwzy;
			Swizzle4D<T, dimension, 1, 0, 2, 3> yxzw;
			Swizzle4D<T, dimension, 1, 0, 3, 2> yxwz;
			Swizzle4D<T, dimension, 1, 2, 0, 3> yzxw;
			Swizzle4D<T, dimension, 1, 2, 3, 0> yzwx;
			Swizzle4D<T, dimension, 1, 3, 0, 2> ywxz;
			Swizzle4D<T, dimension, 1, 3, 2, 0> ywzx;
			Swizzle4D<T, dimension, 2, 0, 1, 3> zxyw;
			Swizzle4D<T, dimension, 2, 0, 3, 1> zxwy;
			Swizzle4D<T, dimension, 2, 1, 0, 3> zyxw;
			Swizzle4D<T, dimension, 2, 1, 3, 0> zywx;
			Swizzle4D<T, dimension, 2, 3, 0, 1> zwxy;
			Swizzle4D<T, dimension, 2, 3, 1, 0> zwyx;
			Swizzle4D<T, dimension, 3, 0, 1, 2> wxyz;
			Swizzle4D<T, dimension, 3, 0, 2, 1> wxzy;
			Swizzle4D<T, dimension, 3, 1, 0, 2> wyxz;
			Swizzle4D<T, dimension, 3, 1, 2, 0> wyzx;
			Swizzle4D<T, dimension, 3, 2, 0, 1> wzxy;
			Swizzle4D<T, dimension, 3, 2, 1, 0> wzyx;

			Swizzle1D<T, dimension, 0> r;
			Swizzle1D<T, dimension, 1> g;
			Swizzle1D<T, dimension, 2> b;
			Swizzle1D<T, dimension, 3> a;
			Swizzle2D<T, dimension, 0, 1> rg;
			Swizzle2D<T, dimension, 1, 0> gr;
			Swizzle2D<T, dimension, 0, 2> rb;
			Swizzle2D<T, dimension, 2, 0> br;
			Swizzle2D<T, dimension, 0, 3> ra;
			Swizzle2D<T, dimension, 3, 0> ar;
			Swizzle2D<T, dimension, 1, 2> gb;
			Swizzle2D<T, dimension, 2, 1> bg;
			Swizzle2D<T, dimension, 1, 3> ga;
			Swizzle2D<T, dimension, 3, 1> ag;
			Swizzle2D<T, dimension, 2, 3> ba;
			Swizzle2D<T, dimension, 3, 2> ab;
			Swizzle3D<T, dimension, 0, 1, 2> rgb;
			Swizzle3D<T, dimension, 0, 2, 1> rbg;
			Swizzle3D<T, dimension, 1, 0, 2> grb;
			Swizzle3D<T, dimension, 1, 2, 0> gbr;
			Swizzle3D<T, dimension, 2, 0, 1> brg;
			Swizzle3D<T, dimension, 2, 1, 0> bgr;
			Swizzle3D<T, dimension, 0, 1, 3> rga;
			Swizzle3D<T, dimension, 0, 3, 1> rag;
			Swizzle3D<T, dimension, 1, 0, 3> gra;
			Swizzle3D<T, dimension, 1, 3, 0> gar;
			Swizzle3D<T, dimension, 3, 0, 1> arg;
			Swizzle3D<T, dimension, 3, 1, 0> agr;
			Swizzle3D<T, dimension, 1, 2, 3> gba;
			Swizzle3D<T, dimension, 1, 3, 2> gab;
			Swizzle3D<T, dimension, 2, 1, 3> bga;
			Swizzle3D<T, dimension, 2, 3, 1> bag;
			Swizzle3D<T, dimension, 3, 1, 2> agb;
			Swizzle3D<T, dimension, 3, 2, 1> abg;
			Swizzle3D<T, dimension, 0, 2, 3> rba;
			Swizzle3D<T, dimension, 0, 3, 2> rab;
			Swizzle3D<T, dimension, 2, 0, 3> bra;
			Swizzle3D<T, dimension, 2, 3, 0> bar;
			Swizzle3D<T, dimension, 3, 0, 2> arb;
			Swizzle3D<T, dimension, 3, 2, 0> abr;
			Swizzle4D<T, dimension, 0, 1, 2, 3> rgba;
			Swizzle4D<T, dimension, 0, 1, 3, 2> rgab;
			Swizzle4D<T, dimension, 0, 2, 1, 3> rbga;
			Swizzle4D<T, dimension, 0, 2, 3, 1> rbag;
			Swizzle4D<T, dimension, 0, 3, 1, 2> ragb;
			Swizzle4D<T, dimension, 0, 3, 2, 1> rabg;
			Swizzle4D<T, dimension, 1, 0, 2, 3> grba;
			Swizzle4D<T, dimension, 1, 0, 3, 2> grab;
			Swizzle4D<T, dimension, 1, 2, 0, 3> gbra;
			Swizzle4D<T, dimension, 1, 2, 3, 0> gbar;
			Swizzle4D<T, dimension, 1, 3, 0, 2> garb;
			Swizzle4D<T, dimension, 1, 3, 2, 0> gabr;
			Swizzle4D<T, dimension, 2, 0, 1, 3> brga;
			Swizzle4D<T, dimension, 2, 0, 3, 1> brag;
			Swizzle4D<T, dimension, 2, 1, 0, 3> bgra;
			Swizzle4D<T, dimension, 2, 1, 3, 0> bgar;
			Swizzle4D<T, dimension, 2, 3, 0, 1> barg;
			Swizzle4D<T, dimension, 2, 3, 1, 0> bagr;
			Swizzle4D<T, dimension, 3, 0, 1, 2> argb;
			Swizzle4D<T, dimension, 3, 0, 2, 1> arbg;
			Swizzle4D<T, dimension, 3, 1, 0, 2> agrb;
			Swizzle4D<T, dimension, 3, 1, 2, 0> agbr;
			Swizzle4D<T, dimension, 3, 2, 0, 1> abrg;
			Swizzle4D<T, dimension, 3, 2, 1, 0> abgr;
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

	// 3D叉积
	template<typename T>
	Vector<T, 3> cross(const Vector<T, 3>& a, const Vector<T, 3>& b)
	{
		return {
			a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x
		};
	}

	// 2D叉积（返回标量）
	template<typename T>
	T cross2D(const Vector<T, 2>& a, const Vector<T, 2>& b)
	{
		return a.x * b.y - a.y * b.x;
	}

	// 距离函数
	template<typename T, size_t D>
	T distance(const Vector<T, D>& a, const Vector<T, D>& b, bool dimensionalityReduction = true)
	{
		return (b - a).length(dimensionalityReduction);
	}

	template<typename T, size_t D>
	T distanceSquared(const Vector<T, D>& a, const Vector<T, D>& b, bool dimensionalityReduction = true)
	{
		return (b - a).lengthSquared(dimensionalityReduction);
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

	// 线性插值
	template<typename T, size_t D>
	Vector<T, D> lerp(const Vector<T, D>& a, const Vector<T, D>& b, T t)
	{
		return a.lerp(b, t);
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

	
/*
 * 这一部分内容为Swizzle的定义部分
 * 
 */
	
	template <typename T, size_t dimension, int idx>
	struct Swizzle1D
	{
		Swizzle1D& operator=(const T& e)
		{
			static_assert(idx < dimension, 
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			elem(idx) = e;
			return *this;
		}

		operator T() const
		{
			static_assert(idx < dimension, 
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			return T(elem(idx));
		}
		
		friend std::ostream& operator<<(std::ostream& out, const Swizzle1D& sw)
		{
			return out << T(sw); 
		}
		
	protected:
		float elem(int i) const
		{
			// 注意，它存在于Union之中，尽管自身仅占1个字节，
			// 但是可以转化为[x,y,z,w]的内存布局的float[4]
			return reinterpret_cast<const T*>(this)[i];
		}

		float& elem(int i)
		{
			return reinterpret_cast<T*>(this)[i];
		}
	};

	template <typename T, size_t dimension, int E0, int E1>
	struct Swizzle2D
	{
		Swizzle2D& operator=(const Vector<T,2>& vec)
		{
			static_assert(E0 < dimension && E1 < dimension, 
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			elem(E0) = vec.x;
			elem(E1) = vec.y;
			return *this;
		}

		operator Vector<T,2>() const
		{
			static_assert(E0 < dimension && E1 < dimension,  
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			return Vector<T,2>{elem(E0),elem(E1)};
		}

		friend std::ostream& operator<<(std::ostream& out, const Swizzle2D& sw)
		{
			return out << Vector<T,2>(sw);
		}

	protected:
		float elem(int i) const
		{
			// 注意，它存在于Union之中，尽管自身仅占1个字节，
			// 但是可以转化为[x,y,z,w]的内存布局的float[4]
			return reinterpret_cast<const T*>(this)[i];
		}

		float& elem(int i)
		{
			return reinterpret_cast<T*>(this)[i];
		}
	};

	template <typename T, size_t dimension, int E0, int E1, int E2>
	struct Swizzle3D
	{
		template<bool Enable = (E0 < dimension && E1 < dimension && E2 < dimension)>
		[[deprecated("Swizzle component access out of bounds! Check this line.")]]
		typename std::enable_if_t<!Enable, Swizzle3D&>
		operator=(const Vector<T,3>& vec)
		{
			static_assert(E0 < dimension && E1 < dimension && E2 < dimension, 
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			elem(E0) = vec.x;
			elem(E1) = vec.y;
			elem(E2) = vec.z;
			return *this;
		}
		
		operator Vector<T,3>() const
		{
			static_assert(E0 < dimension && E1 < dimension && E2 < dimension, 
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			return Vector<T,3>{elem(E0),elem(E1), elem(E2)};
		}
		
		friend std::ostream& operator<<(std::ostream& out, const Swizzle3D& sw)
		{
			static_assert(E0 < dimension && E1 < dimension && E2 < dimension, 
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			return out << Vector<T,3>(sw);
		}

	protected:
		float elem(int i) const
		{
			// 注意，它存在于Union之中，尽管自身仅占1个字节，
			// 但是可以转化为[x,y,z,w]的内存布局的float[4]
			return reinterpret_cast<const T*>(this)[i];
		}

		float& elem(int i)
		{
			return reinterpret_cast<T*>(this)[i];
		}
	};

	template <typename T, size_t dimension, int E0, int E1, int E2, int E3>
	struct Swizzle4D
	{
		Swizzle4D& operator=(const Vector<T,4>& vec)
		{
			static_assert(E0 < dimension && E1 < dimension && E2 < dimension && E3 < dimension,
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			elem(E0) = vec.x;
			elem(E1) = vec.y;
			elem(E2) = vec.z;
			elem(E3) = vec.w;
			return *this;
		}

		operator Vector<T,4>() const
		{
			static_assert(E0 < dimension && E1 < dimension && E2 < dimension && E3 < dimension, 
			"Swizzle index out of bounds! "
			"You're trying to access a component that doesn't exist in this vector.");
			return Vector<T,4>{elem(E0),elem(E1), elem(E2), elem(E3)};
		}

		friend std::ostream& operator<<(std::ostream& out, const Swizzle4D& sw)
		{
			return out << Vector<T,4>(sw);
		}

	protected:
		float elem(int i) const
		{
			// 注意，它存在于Union之中，尽管自身仅占1个字节，
			// 但是可以转化为[x,y,z,w]的内存布局的float[4]
			return reinterpret_cast<const T*>(this)[i];
		}

		float& elem(int i)
		{
			return reinterpret_cast<T*>(this)[i];
		}
	};


/*
 * 常用的向量类型定义
 * 
 */

	using Vector2 = Vector<float, 2>;
	using Vector3 = Vector<float, 3>;
	using Vector4 = Vector<float, 4>;
	
	using Vector2f = Vector<float, 2>;
	using Vector3f = Vector<float, 3>;
	using Vector4f = Vector<float, 4>;

	using Vector2d = Vector<double, 2>;
	using Vector3d = Vector<double, 3>;
	using Vector4d = Vector<double, 4>;

	using Vector2i = Vector<int, 2>;
	using Vector3i = Vector<int, 3>;
	using Vector4i = Vector<int, 4>;

}
