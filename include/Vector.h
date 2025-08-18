#pragma once

#include "EasyMathAPI.h"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <array>
#include <cassert>

#include "Common.h"

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

	template<typename T, size_t dimension>
	class Vector
	{
	public:
		static_assert(std::is_same_v<T, float> || std::is_same_v<T, double> || std::is_same_v<T, long double>, "Member of vector must be float, double or long double.");

		static_assert(dimension > 0, "vector dimension must be positive.");

		Vector()
			:data({})
		{
		};

		Vector(const T& num)
		{
			for (auto it = data.begin(); it != data.end(); ++it)
			{
				*it = num;
			}
		};

		Vector(std::initializer_list<T> InitializeList)
		{
			assert(InitializeList.size() == dimension && "[Vector] Invalid initializer list size");

			std::copy(InitializeList.begin(), InitializeList.end(), data.begin());
		}

		T& operator[](size_t idx)
		{
			limit_max(idx, dimension - 1);
			return data[idx];
		}

		const T& operator[](size_t idx) const
		{
			limit_max(idx, dimension - 1);
			return data[idx];
		}

		T& operator[](VectorFilterDimension1D d)
		{
			uint8_t idx = (uint8_t)d % 4;

			assert(idx < dimension && "dimension index attached the vector border.");

			return data[(size_t)idx];
		}

		const T& operator[](VectorFilterDimension1D d) const
		{
			uint8_t idx = (uint8_t)d % 4;

			assert(idx < dimension && "dimension index attached the vector border.");

			return data[(size_t)idx];
		}

		Vector<T, 2> operator[](VectorFilterDimension2D d)
		{
			assert(dimension > 1 && "vector dimension doesn't support double swizzle.");

			uint8_t idx = (uint8_t)d % 6;

			if (idx < 1)
			{
				return Vector<T, 2>{data[0], data[1]};
			}

			if (idx < 3)
			{
				assert(dimension > 2 && "dimension index attached the vector border.");

				if (idx == 1)
				{
					return Vector<T, 2>{data[0], data[2]};
				}

				if (idx == 2)
				{
					return Vector<T, 2>{data[1], data[2]};
				}
			}


			assert(dimension > 3 && "dimension index attached the vector border.");

			if (idx == 3)
			{
				return Vector<T, 2>{data[0], data[3]};
			}

			if (idx == 4)
			{
				return Vector<T, 2>{data[1], data[3]};
			}

			return Vector<T, 2>{data[2], data[3]};

		}

		const Vector<T, 2> operator[](VectorFilterDimension2D d) const
		{
			assert(dimension > 1 && "vector dimension doesn't support double swizzle.");

			uint8_t idx = (uint8_t)d % 6;

			if (idx < 1)
			{
				return Vector<T, 2>{data[0], data[1]};
			}

			if (idx < 3)
			{
				assert(dimension > 2 && "dimension index attached the vector border.");

				if (idx == 1)
				{
					return Vector<T, 2>{data[0], data[2]};
				}

				if (idx == 2)
				{
					return Vector<T, 2>{data[1], data[2]};
				}
			}

			assert(dimension > 3 && "dimension index attached the vector border.");

			if (idx == 3)
			{
				return Vector<T, 2>{data[0], data[3]};
			}

			if (idx == 4)
			{
				return Vector<T, 2>{data[1], data[3]};
			}

			return Vector<T, 2>{data[2], data[3]};

		}

		Vector<T, 3> operator[](VectorFilterDimension3D d)
		{
			assert(dimension > 2 && "vector dimension doesn't support triple swizzle.");

			uint8_t idx = (uint8_t)d % 3;

			if (idx < 1)
			{
				return Vector<T, 3>{data[0], data[1], data[2]};
			}

			assert(dimension > 3 && "dimension index attached the vector border.");

			if (idx == 1)
			{
				return Vector<T, 3>{data[0], data[1], data[3]};
			}

			return Vector<T, 3>{data[1], data[2], data[3]};
		}

		const Vector<T, 3> operator[](VectorFilterDimension3D d) const
		{
			assert(dimension > 2 && "vector dimension doesn't support triple swizzle.");

			uint8_t idx = (uint8_t)d % 3;

			if (idx < 1)
			{
				return Vector<T, 3>{data[0], data[1], data[2]};
			}

			assert(dimension > 3 && "dimension index attached the vector border.");

			if (idx == 1)
			{
				return Vector<T, 3>{data[0], data[1], data[3]};
			}

			return Vector<T, 3>{data[1], data[2], data[3]};
		}

		T Length(bool dimensionalityReduction = true) const
		{
			T sum = T(0);

			for (size_t idx = 0; idx < dimension; ++idx)
			{
				if (dimensionalityReduction && idx > 2)
				{
					break;
				}

				sum += data[idx] * data[idx];
			}

			return std::sqrt(sum);
		}

		Vector<T, dimension> Normalize(bool dimensionalityReduction = true) const
		{
			T len = Length(dimensionalityReduction);

			Vector<T, dimension> result;

			for (size_t idx = 0; idx < dimension; ++idx)
			{
				if (dimensionalityReduction && idx > 2)
				{
					result[idx] = data[idx];
					continue;
				}

				if (len == 0)
				{
					result[idx] = 0;
					continue;
				}

				result[idx] = data[idx] / len;
			}

			return result;
		}

		static constexpr size_t Dimension() { return dimension; }

		Vector<T, dimension>& operator=(const Vector<T, dimension>& other)
		{
			for (size_t idx = 0; idx < dimension; ++idx)
			{
				this->data[idx] = other[idx];
			}

			return *this;
		}

	protected:
		std::array<T, dimension> data;

	private:

	};



	template<typename T, size_t dimension_var>
	Vector<T, 3> Cross(const Vector<T, dimension_var>& vectorA, const Vector<T, dimension_var>& vectorB)
	{
		static_assert(dimension_var == 2 || dimension_var == 3, "Cross vector must be Vector2 or Vector3.");

		if (dimension_var == 2)
		{
			return Vector<T, 3>{ 0, 0, vectorA[x] * vectorB[y] - vectorA[y] * vectorB[x] };
		}

		return Vector<T, 3>{ vectorA[y] * vectorB[z] - vectorA[z] * vectorB[y],
			vectorA[z] * vectorB[x] - vectorA[x] * vectorB[z],
			vectorA[x] * vectorB[y] - vectorA[y] * vectorB[x] };
	}

	template<typename T, size_t dimension_var>
	T Dot(const Vector<T, dimension_var>& vectorA, const Vector<T, dimension_var>& vectorB, bool dimensionalityReduction = true)
	{
		T result = 0;

		for (size_t idx = 0; idx < dimension_var; ++idx)
		{
			if (dimensionalityReduction && idx > 2)
			{
				break;
			}

			result += vectorA[idx] * vectorB[idx];
		}

		return result;
	}

	template<typename T, size_t D>
	std::ostream& operator<<(std::ostream& out, const Vector<T, D>& vector)
	{
		out << "(";

		for (size_t idx = 0; idx < D; ++idx)
		{
			out << vector[idx];

			if (idx != D - 1)
			{
				out << ",";
			}
		}

		out << ")";

		return out;
	}

	template<typename T, size_t D>
	Vector<T, D>& operator+=(Vector<T, D>& vector, T scalar)
	{
		for (size_t idx = 0; idx < D; ++idx)
		{
			vector[idx] += scalar;
		}

		return vector;
	}

	template<typename T, size_t D>
	Vector<T, D>& operator+=(Vector<T, D>& vector, const Vector<T, D>& other)
	{
		for (size_t idx = 0; idx < D; ++idx)
		{
			vector[idx] += other[idx];
		}

		return vector;
	}

	template<typename T, size_t D>
	Vector<T, D>& operator*=(Vector<T, D>& vector, T scalar)
	{
		for (size_t idx = 0; idx < D; ++idx)
		{
			vector[idx] *= scalar;
		}

		return vector;
	}

	template<typename T, size_t D>
	Vector<T, D>& operator*=(Vector<T, D>& vector, const Vector<T, D>& other)
	{
		for (size_t idx = 0; idx < D; ++idx)
		{
			vector[idx] *= other[idx];
		}

		return vector;
	}

	template<typename T, size_t D>
	Vector<T, D>& operator-=(Vector<T, D>& vector, T scalar)
	{
		for (size_t idx = 0; idx < D; ++idx)
		{
			vector[idx] -= scalar;
		}

		return vector;
	}

	template<typename T, size_t D>
	Vector<T, D>& operator-=(Vector<T, D>& vector, const Vector<T, D>& other)
	{
		for (size_t idx = 0; idx < D; ++idx)
		{
			vector[idx] -= other[idx];
		}

		return vector;
	}

	template<typename T, size_t D>
	Vector<T, D>& operator/=(Vector<T, D>& vector, T scalar)
	{
		for (size_t idx = 0; idx < D; ++idx)
		{
			vector[idx] /= scalar;
		}

		return vector;
	}

	template<typename T, size_t D>
	Vector<T, D>& operator/=(Vector<T, D>& vector, const Vector<T, D>& other)
	{
		for (size_t idx = 0; idx < D; ++idx)
		{
			vector[idx] /= other[idx];
		}

		return vector;
	}

	template<typename T, size_t D>
	Vector<T, D> operator*(const Vector<T, D>& vectorA, const Vector<T, D>& vectorB)
	{
		Vector<T, D> result;

		for (size_t idx = 0; idx < D; ++idx)
		{
			result[idx] = vectorA[idx] * vectorB[idx];
		}

		return result;
	}

	template<typename T, size_t D>
	Vector<T, D> operator*(const Vector<T, D>& vector, T scalar)
	{
		Vector<T, D> result;

		for (size_t idx = 0; idx < D; ++idx)
		{
			result[idx] = scalar * vector[idx];
		}

		return result;
	}

	template<typename T, size_t D>
	Vector<T, D> operator+(const Vector<T, D>& vectorA, const Vector<T, D>& vectorB)
	{
		Vector<T, D> result;

		for (size_t idx = 0; idx < D; ++idx)
		{
			result[idx] = vectorA[idx] + vectorB[idx];
		}

		return result;
	}

	template<typename T, size_t D>
	Vector<T, D> operator+(const Vector<T, D>& vector, T scalar)
	{
		Vector<T, D> result;

		for (size_t idx = 0; idx < D; ++idx)
		{
			result[idx] = scalar + vector[idx];
		}

		return result;
	}

	template<typename T, size_t D>
	Vector<T, D> operator-(const Vector<T, D>& vectorA, const Vector<T, D>& vectorB)
	{
		Vector<T, D> result;

		for (size_t idx = 0; idx < D; ++idx)
		{
			result[idx] = vectorA[idx] - vectorB[idx];
		}

		return result;
	}

	template<typename T, size_t D>
	Vector<T, D> operator-(const Vector<T, D>& vector, T scalar)
	{
		Vector<T, D> result;

		for (size_t idx = 0; idx < D; ++idx)
		{
			result[idx] = vector[idx] - scalar;
		}

		return result;
	}

	template<typename T, size_t D>
	Vector<T, D> operator/(const Vector<T, D>& vectorA, const Vector<T, D>& vectorB)
	{
		Vector<T, D> result;

		for (size_t idx = 0; idx < D; ++idx)
		{
			result[idx] = vectorA[idx] / vectorB[idx];
		}

		return result;
	}

	template<typename T, size_t D>
	Vector<T, D> operator/(const Vector<T, D>& vector, T scalar)
	{
		Vector<T, D> result;

		for (size_t idx = 0; idx < D; ++idx)
		{
			result[idx] = vector[idx] / scalar;
		}

		return result;
	}

	typedef Vector<float, 2> Vector2;
	typedef Vector<float, 3> Vector3;
	typedef Vector<float, 4> Vector4;

	EASYMATH_API bool NearZero(const Vector2& src);
	EASYMATH_API bool NearZero(const Vector3& src);

	EASYMATH_API Vector2 DegreeToRadians(const Vector2& degree);
	EASYMATH_API Vector3 DegreeToRadians(const Vector3& degree);

	EASYMATH_API Vector2 RadiansToDegree(Vector2& radians);
	EASYMATH_API Vector3 RadiansToDegree(Vector3& radians);

}