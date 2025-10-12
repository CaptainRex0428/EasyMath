#pragma once

#include "EasyMathAPI.h"
#include "Vector.h"

namespace EM
{
    template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    class Quaternion : public Vector<T,4>
    {
    public:
        Quaternion():Vector<T,4>({0,0,0,0}){}
        
        Quaternion(T x, T y, T z, T w):Vector<T,4>({x,y,z,w}){}
        Quaternion(T w, Vector<T, 3> xyz):Vector<T,4>({xyz[x],xyz[y],xyz[z], w}){}
        Quaternion(Vector<T, 3> xyz, T w):Vector<T,4>({xyz[x],xyz[y],xyz[z], w}){}

        // 输出
        virtual std::ostream& print(std::ostream& out) const override
        {
            out << "[" << (*this)[w] << "," << (*this)[x] << "," << (*this)[y] << "," << (*this)[z] << "]";
            return out;
        }
    
    };

    // 类型别名
    using QuaternionF = Quaternion<float>;
    using QuaternionD = Quaternion<double>;
}


