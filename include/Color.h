#pragma once

#include "EasyMathAPI.h"
#include "EMConst.h"
#include "Vector.h"


namespace EM
{
    enum class LuminanceStandard
    {
        Rec601,    // 0.299, 0.587, 0.114 (NTSC)
        Rec709,    // 0.2126, 0.7152, 0.0722 (HDTV)
        Rec2020    // 0.2627, 0.6780, 0.0593 (UHD)
    };
    
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
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

        static Vector<T, 4> getRGBLuminanceContribute(LuminanceStandard ls =  LuminanceStandard::Rec709)
        {
            Vector<T, 4> result;
            
            switch (ls)
            {
            case LuminanceStandard::Rec2020:
                result = Vector<T, 4>{0.2627f, 0.6780f, 0.0593f, 1.f};
                break;
            case LuminanceStandard::Rec709:
                 result = Vector<T, 4>{0.2126f, 0.7152f, 0.0722f, 1.f};
                break;
            case LuminanceStandard::Rec601:
            default:
                result = Vector<T, 4>{0.299f, 0.587f, 0.114f, 1.f};
                break;
            }

            return result;
        }
        
        T luminance(LuminanceStandard ls = LuminanceStandard::Rec709)
        {
            return dot((*this).rgb,getRGBLuminanceContribute(ls).rgb);
        }
    };
}

