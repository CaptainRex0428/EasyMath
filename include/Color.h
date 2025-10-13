#pragma once

#include "EasyMathAPI.h"
#include "EMConst.h"
#include "Vector.h"


namespace EM
{
    template <typename T,typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    class Color : public Vector<T, 4>
    {
    public:
        Color()
            :Vector<T, 4>()
        {
        }

        explicit Color(const T& value):Vector<T, 4>(value)
        {
        }

        Color(std::initializer_list<T> InitializeList):Vector<T, 4>(InitializeList)
        {
        }

        Color(const Color&) = default;

        Color(Color&&) = default;

        virtual ~Color() = default;
        
        // 颜色通道对亮度的贡献权重
        static Vector<T,4> ColorLuminanceContribute()
        {
            return Vector<T, 4>({0.299f, 0.587f, 0.114f, 0});
        }

        T luminance()
        {
            return dot(*this,ColorLuminanceContribute());
        }
    };
}

