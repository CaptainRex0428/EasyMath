#pragma once

#include "EasyMathAPI.h"
#include "Common.h"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <array>
#include <cassert>
#include <algorithm>

namespace EM
{
	enum VectorFilterDimension1D : uint8_t
	{
		x = 0, y, z, w,
		r, g, b, a,
		X, Y, Z, W,
		R, G, B, A
	};

	enum VectorFilterDimension2D : uint8_t
	{
		xy = 0, xz, yz, xw, yw, zw,
		rg, rb, gb, ra, ga, ba,
		XY, XZ, YZ, XW, YW, ZW,
		RG, RB, GB, RA, GA, BA
	};

	enum VectorFilterDimension3D : uint8_t
	{
		xyz = 0, xyw, yzw,
		rgb, rga, gba,
		XYZ, XYW, YZW,
		RGB, RGA, GBA
	};

	template <typename T, size_t rows, size_t cols, typename>
	class Matrix;

	template<typename T, typename>
	class Quaternion;

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

		Vector(Quaternion<T, std::enable_if_t<std::is_arithmetic_v<T>>>&& Q)
		{
			static_assert(dimension == 4, "dimension must be 4 as a quaternion.");

			this->data[x] = Q.swizzle(x);
			this->data[y] = Q.swizzle(y);
			this->data[z] = Q.swizzle(z);
			this->data[w] = Q.swizzle(w);
		}

		virtual ~Vector() = default;

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

		// Swizzle
		virtual const T& operator[](VectorFilterDimension1D d) const
		{
			uint8_t idx = (uint8_t)d % 4;
			assert(idx < dimension && "Swizzle index out of bounds");
			return data[idx];
		}
		
		virtual T& operator[](VectorFilterDimension1D d)
		{
			uint8_t idx = (uint8_t)d % 4;
			assert(idx < dimension && "Swizzle index out of bounds");
			return data[idx];
		}
		
		const T& swizzle(VectorFilterDimension1D d) const
		{
			uint8_t idx = (uint8_t)d % 4;
			assert(idx < dimension && "Swizzle index out of bounds");
			return data[idx];
		}

		T& swizzle(VectorFilterDimension1D d)
		{
			uint8_t idx = (uint8_t)d % 4;
			assert(idx < dimension && "Swizzle index out of bounds");
			return data[idx];
		}

		// 2D Swizzle
		Vector<T, 2> swizzle(VectorFilterDimension2D d) const
		{
			assert(dimension >= 2 && "Vector dimension too small for 2D swizzle");
			uint8_t idx = (uint8_t)d % 6;

			switch (idx) 
			{
				case 0: return { data[0], data[1] };  // xy
				case 1: assert(dimension >= 3); return { data[0], data[2] };  // xz
				case 2: assert(dimension >= 3); return { data[1], data[2] };  // yz
				case 3: assert(dimension >= 4); return { data[0], data[3] };  // xw
				case 4: assert(dimension >= 4); return { data[1], data[3] };  // yw
				case 5: assert(dimension >= 4); return { data[2], data[3] };  // zw
				default: return { data[0], data[1] };
			}
		}

		// 3D Swizzle
		Vector<T, 3> swizzle(VectorFilterDimension3D d) const
		{
			assert(dimension >= 3 && "Vector dimension too small for 3D swizzle");
			uint8_t idx = static_cast<uint8_t>(d) % 3;

			switch (idx) {
			case 0: return { data[0], data[1], data[2] };  // xyz
			case 1: assert(dimension >= 4); return { data[0], data[1], data[3] };  // xyw
			case 2: assert(dimension >= 4); return { data[1], data[2], data[3] };  // yzw
			default: return { data[0], data[1], data[2] };
			}
		}

		// 转换为矩阵
		virtual Matrix<T, dimension, 1, std::enable_if_t<std::is_arithmetic_v<T>>> toColMatrix()
		{
			Matrix<T, dimension, 1> result{};
			for (size_t i = 0; i < dimension; ++i) 
			{
				result(i, 0) = this->data[i];
			}
			return result;
		}

		virtual Matrix<T, 1, dimension, std::enable_if_t<std::is_arithmetic_v<T>>> toRowMatrix()
		{
			Matrix<T, 1, dimension> result{};
			for (size_t i = 0; i < dimension; ++i) 
			{
				result(0, i) = this->data[i];
			}
			return result;
		}

		// 成员函数：复合赋值运算符
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

		// 向量运算成员函数

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

		// 图形渲染专用函数
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

		virtual bool isNormalized(T epsilon = T{ NEARZERO_THRESHOLD }, bool dimensionalityReduction = true) const noexcept
		{
			return std::abs(length(dimensionalityReduction) - T{ 1 }) <= epsilon;
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
			Matrix<T, 4, 4> result = MTXIdentity<T, 4>();
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

		// 数据访问
		virtual T* Data() noexcept { return data.data(); }
		virtual const T* Data() const noexcept { return data.data(); }

		virtual T* begin() noexcept { return data.data(); }
		virtual const T* begin() const noexcept { return data.data(); }

		virtual T* end() noexcept { return data.data() + dimension; }
		virtual const T* end() const noexcept { return data.data() + dimension; }

		// 基础信息
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

		// 输出
		virtual std::ostream& print(std::ostream& out) const
		{
			out << "(";
			for (size_t i = 0; i < dimension; ++i) 
			{
				out << (*this).data[i];
				if (i < dimension - 1) {
					out << ", ";
				}
			}
			out << ") (V mode: xyzw)";
			return out;
		}

	private:
		// 友元函数：二元运算符
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

	protected:
		std::array<T, dimension> data;

	};

	template<typename T, size_t dimension>
	std::ostream& operator<<(std::ostream& out, const Vector<T,dimension>& vec)
	{
		return vec.print(out);
	}

	// 全局向量函数
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
			a.swizzle(y) * b.swizzle(z) - a.swizzle(z) * b.swizzle(y),
			a.swizzle(z) * b.swizzle(x) - a.swizzle(x) * b.swizzle(z),
			a.swizzle(x) * b.swizzle(y) - a.swizzle(y) * b.swizzle(x)
		};
	}

	// 2D叉积（返回标量）
	template<typename T>
	T cross2D(const Vector<T, 2>& a, const Vector<T, 2>& b)
	{
		return a.swizzle(x) * b.swizzle(y) - a.swizzle(y) * b.swizzle(x);
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

	// 常用的向量类型定义
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
