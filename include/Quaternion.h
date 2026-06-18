#pragma once

#include "EasyMathAPI.h"
#include "Vector.h"

namespace EasyMath
{
    template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    class Quaternion : public Vector<T,4>
    {
    public:
        Quaternion():Vector<T,4>({0,0,0,0}){}
        
        Quaternion(T x, T y, T z, T w):Vector<T,4>({x,y,z,w}){}
        Quaternion(T w, Vector<T, 3> xyz):Vector<T,4>({xyz.x,xyz.y,xyz.z, w}){}
        Quaternion(Vector<T, 3> xyz, T w):Vector<T,4>({xyz.x,xyz.y,xyz.z, w}){}
        Quaternion(Vector<T,4> xyzw): Vector<T,4>(xyzw){}
        
        virtual T length(bool dimensionalityReduction = false) const noexcept override
        {
            assert(!dimensionalityReduction && "Quaternion does not support dimensionalityReduction");
            
            return std::sqrt(lengthSquared(dimensionalityReduction));
        }

        virtual constexpr T lengthSquared(bool dimensionalityReduction = false) const noexcept override
        {
            assert(!dimensionalityReduction && "Quaternion does not support dimensionalityReduction");
            
            T sum = T(0);
            for (size_t idx = 0; idx < 4; ++idx)
            {
                sum += this->data[idx] * this->data[idx];
            }
            return sum;
        }

        virtual void normalize(bool dimensionalityReduction = false) override
        {
            assert(!dimensionalityReduction && "Quaternion does not support dimensionalityReduction");
            *this = normalized();
        }

        virtual Vector<T, 4> normalized(bool) const override
        {
            return normalized();
        }
        
        virtual Quaternion normalized() const
        {
            T len = this->length();
            if (len == T{ 0 }) return Quaternion();
        
            Quaternion result;
            for (size_t idx = 0; idx < 4; ++idx)
            {
                result[idx] = (*this)[idx] / len;
            }
            return result;
        }
        
        // 输出
        virtual std::ostream& print(std::ostream& out) const override
        {
            out << "[" << this->w << "," << this->x << "," << this->y << "," << this->z << "] (Q mode: wxyz)";
            return out;
        }
    
    };

    // 类型别名
    using QuaternionF = Quaternion<float>;
    using QuaternionD = Quaternion<double>;
}


