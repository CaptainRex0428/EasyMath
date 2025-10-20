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

    template <typename T, bool bsRGB ,bool bRGBMode, typename>
    class Color;

    template <typename T>
    using sRGBColor = Color<T, true, true,  std::enable_if_t<std::is_arithmetic_v<T>>>;

    template <typename T>
    using LinearColor = Color<T, true, false, std::enable_if_t<std::is_arithmetic_v<T>>>;

    template <typename T>
    using HSV = Color<T, false, false,  std::enable_if_t<std::is_arithmetic_v<T>>>;

    template <typename T>
    using HSL = Color<T, false, false,  std::enable_if_t<std::is_arithmetic_v<T>>>;

    template <typename T>
    using HSI = Color<T, false, false,  std::enable_if_t<std::is_arithmetic_v<T>>>;
    
    template <typename T,bool isRGBMode = true, bool bsRGB = false , typename = std::enable_if_t<std::is_arithmetic_v<T>>>
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

        T luminance(LuminanceStandard ls = LuminanceStandard::Rec709)
        {
            return dot((*this).rgb,getRGBLuminanceContribute(ls).rgb);
        }

        /**
         * 将sRGB颜色转换为线性颜色空间
         * @return 线性颜色空间的Color对象
         */
        template <bool rgbmode = isRGBMode, bool srgb = bsRGB>
        typename std::enable_if_t<rgbmode && bsRGB, LinearColor<T>>
        toLinear() const
        {
            Color result;
            for (size_t i = 0; i < 3; ++i)
            {
                T value = (*this)[i];
                if (value <= T(0.04045))
                {
                    result[i] = value / T(12.92);
                }
                else
                {
                    result[i] = std::pow((value + T(0.055)) / T(1.055), T(2.4));
                }
            }
            result[3] = (*this)[3]; // Alpha通道保持不变
            return result;
        }

        /**
         * 将线性颜色空间转换为sRGB
         * @return sRGB颜色空间的Color对象
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode && (!bsRGB), sRGBColor<T>>
        toSRGB() const
        {
            Color result;
            for (size_t i = 0; i < 3; ++i)
            {
                T value = (*this)[i];
                if (value <= T(0.0031308))
                {
                    result[i] = value * T(12.92);
                }
                else
                {
                    result[i] = T(1.055) * std::pow(value, T(1.0 / 2.4)) - T(0.055);
                }
            }
            result[3] = (*this)[3]; // Alpha通道保持不变
            return result;
        }

        /**
         * 将RGB颜色转换为HSV颜色空间
         * @return HSV格式的Vector3 (H: [0,360], S: [0,1], V: [0,1])
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode, HSV<T>>
        toHSV() const
        {
            T r = (*this)[0];
            T g = (*this)[1];
            T b = (*this)[2];

            T max = std::max({r, g, b});
            T min = std::min({r, g, b});
            T delta = max - min;

            Vector<T, 3> hsv;
            
            // Hue计算
            if (delta < T(1e-6))
            {
                hsv[0] = T(0);
            }
            else if (max == r)
            {
                hsv[0] = T(60) * (std::fmod((g - b) / delta, T(6)));
            }
            else if (max == g)
            {
                hsv[0] = T(60) * ((b - r) / delta + T(2));
            }
            else
            {
                hsv[0] = T(60) * ((r - g) / delta + T(4));
            }
            
            if (hsv[0] < T(0))
            {
                hsv[0] += T(360);
            }

            // Saturation计算
            hsv[1] = (max < T(1e-6)) ? T(0) : (delta / max);

            // Value计算
            hsv[2] = max;

            return hsv;
        }

        /**
         * 将RGB颜色转换为HSL颜色空间
         * @return HSL格式的Vector3 (H: [0,360], S: [0,1], L: [0,1])
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode, HSL<T>>
        toHSL() const
        {
            T r = (*this)[0];
            T g = (*this)[1];
            T b = (*this)[2];

            T max = std::max({r, g, b});
            T min = std::min({r, g, b});
            T delta = max - min;

            Vector<T, 3> hsl;
            
            // Lightness计算
            hsl[2] = (max + min) / T(2);

            // Saturation计算
            if (delta < T(1e-6))
            {
                hsl[0] = T(0);
                hsl[1] = T(0);
            }
            else
            {
                // Hue计算（与HSV相同）
                if (max == r)
                {
                    hsl[0] = T(60) * (std::fmod((g - b) / delta, T(6)));
                }
                else if (max == g)
                {
                    hsl[0] = T(60) * ((b - r) / delta + T(2));
                }
                else
                {
                    hsl[0] = T(60) * ((r - g) / delta + T(4));
                }
                
                if (hsl[0] < T(0))
                {
                    hsl[0] += T(360);
                }

                // Saturation (HSL定义)
                hsl[1] = delta / (T(1) - std::abs(T(2) * hsl[2] - T(1)));
            }

            return hsl;
        }

        /**
         * 将RGB颜色转换为HSI颜色空间
         * @return HSI格式的Vector3 (H: [0,360], S: [0,1], I: [0,1])
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode, HSI<T>>
        toHSI() const
        {
            T r = (*this)[0];
            T g = (*this)[1];
            T b = (*this)[2];

            Vector<T, 3> hsi;
            
            // Intensity计算
            hsi[2] = (r + g + b) / T(3);

            // Saturation计算
            T min_rgb = std::min({r, g, b});
            if (hsi[2] < T(1e-6))
            {
                hsi[1] = T(0);
                hsi[0] = T(0);
            }
            else
            {
                hsi[1] = T(1) - min_rgb / hsi[2];
                
                // Hue计算
                if (hsi[1] < T(1e-6))
                {
                    hsi[0] = T(0);
                }
                else
                {
                    T numerator = T(0.5) * ((r - g) + (r - b));
                    T denominator = std::sqrt((r - g) * (r - g) + (r - b) * (g - b));
                    
                    if (denominator < T(1e-6))
                    {
                        hsi[0] = T(0);
                    }
                    else
                    {
                        T theta = std::acos(numerator / denominator);
                        hsi[0] = (b <= g) ? theta * T(180) / T(PI) : T(360) - theta * T(180) / T(PI);
                    }
                }
            }

            return hsi;
        }

    public:

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

        /**
         * 从HSV颜色空间创建RGB颜色
         * @param h 色相 [0, 360]
         * @param s 饱和度 [0, 1]
         * @param v 明度 [0, 1]
         * @param a Alpha通道值 [0, 1]
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode, bool srgb = bsRGB>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, bsRGB>>
        fromHSV(T h, T s, T v, T a = T(1))
        {
            // 将色相规范化到[0, 360)
            h = std::fmod(h, T(360));
            if (h < T(0)) h += T(360);

            T c = v * s;
            T x = c * (T(1) - std::abs(std::fmod(h / T(60), T(2)) - T(1)));
            T m = v - c;

            T r, g, b;

            if (h < T(60))
            {
                r = c; g = x; b = T(0);
            }
            else if (h < T(120))
            {
                r = x; g = c; b = T(0);
            }
            else if (h < T(180))
            {
                r = T(0); g = c; b = x;
            }
            else if (h < T(240))
            {
                r = T(0); g = x; b = c;
            }
            else if (h < T(300))
            {
                r = x; g = T(0); b = c;
            }
            else
            {
                r = c; g = T(0); b = x;
            }

            return Color{r + m, g + m, b + m, a};
        }

        /**
         * 从HSV Vector创建RGB颜色
         * @param hsv HSV格式的Vector3 (H: [0,360], S: [0,1], V: [0,1])
         * @param a Alpha通道值 [0, 1]
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode, bool srgb = bsRGB>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, bsRGB>>
        fromHSV(const Vector<T, 3>& hsv, T a = T(1))
        {
            return fromHSV(hsv[0], hsv[1], hsv[2], a);
        }

        /**
         * 从HSL颜色空间创建RGB颜色
         * @param h 色相 [0, 360]
         * @param s 饱和度 [0, 1]
         * @param l 亮度 [0, 1]
         * @param a Alpha通道值 [0, 1]
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode, bool srgb = bsRGB>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, bsRGB>>
        fromHSL(T h, T s, T l, T a = T(1))
        {
            // 将色相规范化到[0, 360)
            h = std::fmod(h, T(360));
            if (h < T(0)) h += T(360);

            T c = (T(1) - std::abs(T(2) * l - T(1))) * s;
            T x = c * (T(1) - std::abs(std::fmod(h / T(60), T(2)) - T(1)));
            T m = l - c / T(2);

            T r, g, b;

            if (h < T(60))
            {
                r = c; g = x; b = T(0);
            }
            else if (h < T(120))
            {
                r = x; g = c; b = T(0);
            }
            else if (h < T(180))
            {
                r = T(0); g = c; b = x;
            }
            else if (h < T(240))
            {
                r = T(0); g = x; b = c;
            }
            else if (h < T(300))
            {
                r = x; g = T(0); b = c;
            }
            else
            {
                r = c; g = T(0); b = x;
            }

            return Color{r + m, g + m, b + m, a};
        }

        /**
         * 从HSL Vector创建RGB颜色
         * @param hsl HSL格式的Vector3 (H: [0,360], S: [0,1], L: [0,1])
         * @param a Alpha通道值 [0, 1]
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode, bool srgb = bsRGB>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, bsRGB>>
        fromHSL(const Vector<T, 3>& hsl, T a = T(1))
        {
            return fromHSL(hsl[0], hsl[1], hsl[2], a);
        }

        /**
         * 从HSI颜色空间创建RGB颜色
         * @param h 色相 [0, 360]
         * @param s 饱和度 [0, 1]
         * @param i 强度 [0, 1]
         * @param a Alpha通道值 [0, 1]
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode, bool srgb = bsRGB>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, bsRGB>>
        fromHSI(T h, T s, T i, T a = T(1))
        {
            // 将色相规范化到[0, 360)
            h = std::fmod(h, T(360));
            if (h < T(0)) h += T(360);

            T r, g, b;

            // 转换为弧度
            T h_rad = h * T(PI) / T(180);

            if (h < T(120))
            {
                b = i * (T(1) - s);
                r = i * (T(1) + s * std::cos(h_rad) / std::cos(T(PI) / T(3) - h_rad));
                g = T(3) * i - (r + b);
            }
            else if (h < T(240))
            {
                h_rad -= T(2) * T(PI) / T(3);
                r = i * (T(1) - s);
                g = i * (T(1) + s * std::cos(h_rad) / std::cos(T(PI) / T(3) - h_rad));
                b = T(3) * i - (r + g);
            }
            else
            {
                h_rad -= T(4) * T(PI) / T(3);
                g = i * (T(1) - s);
                b = i * (T(1) + s * std::cos(h_rad) / std::cos(T(PI) / T(3) - h_rad));
                r = T(3) * i - (g + b);
            }

            // 钳制到[0, 1]范围
            r = std::max(T(0), std::min(T(1), r));
            g = std::max(T(0), std::min(T(1), g));
            b = std::max(T(0), std::min(T(1), b));

            return Color{r, g, b, a};
        }

        /**
         * 从HSI Vector创建RGB颜色
         * @param hsi HSI格式的Vector3 (H: [0,360], S: [0,1], I: [0,1])
         * @param a Alpha通道值 [0, 1]
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode, bool srgb = bsRGB>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, bsRGB>>
        fromHSI(const Vector<T, 3>& hsi, T a = T(1))
        {
            return fromHSI(hsi[0], hsi[1], hsi[2], a);
        }
    };
    
}

