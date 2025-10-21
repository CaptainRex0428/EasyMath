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

    template <typename T, bool bsRGB ,size_t modeIndex, typename>
    class Color;

    template <typename T>
    using sRGBColor = Color<T, true, 1,  std::enable_if_t<std::is_arithmetic_v<T>>>;

    template <typename T>
    using LinearColor = Color<T, true, 0, std::enable_if_t<std::is_arithmetic_v<T>>>;

    template <typename T>
    using HSV = Color<T, false, 0,  std::enable_if_t<std::is_arithmetic_v<T>>>;

    template <typename T>
    using HSL = Color<T, false, 1,  std::enable_if_t<std::is_arithmetic_v<T>>>;

    template <typename T>
    using HSI = Color<T, false, 2,  std::enable_if_t<std::is_arithmetic_v<T>>>;
    
    template <typename T,bool isRGBMode = true, size_t modeIndex = 0 , typename = std::enable_if_t<std::is_arithmetic_v<T>>>
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

        /**
         * 计算物理正确的亮度（Luminance）
         * 对于sRGB颜色，会先转换到线性空间再计算
         * 对于Linear颜色，直接计算
         * @param ls 亮度标准（Rec601/Rec709/Rec2020）
         * @return 亮度值 [0, 1]
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode, T>
        luminance(LuminanceStandard ls = LuminanceStandard::Rec709)
        {
            if constexpr (modeIndex == 1) // sRGB
            {
                // sRGB需要先转换到线性空间
                auto linear = toLinear();
                return dot(linear.rgb, getRGBLuminanceContribute(ls).rgb);
            }
            
            // 线性空间直接计算
            return dot((*this).rgb, getRGBLuminanceContribute(ls).rgb);
        }

        /**
         * 计算感知亮度（Perceived Luminance）
         * 直接在当前颜色空间（sRGB或Linear）计算
         * 适用于UI显示、文本对比度计算等场景
         * @param ls 亮度标准
         * @return 感知亮度值
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode, T>
        perceivedLuminance(LuminanceStandard ls = LuminanceStandard::Rec709) const
        {
            // 直接在当前颜色空间计算（不转换）
            return dot((*this).rgb, getRGBLuminanceContribute(ls).rgb);
        }
        
        /**
        * 获取相对亮度（Relative Luminance）- W3C WCAG标准
        * 用于计算对比度，符合无障碍访问标准
        * 始终在线性空间计算
        * @return 相对亮度 [0, 1]
        */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode, T>
        relativeLuminance() const
        {
            // WCAG使用Rec.709系数，并要求在线性空间计算
            if constexpr (modeIndex == 1) // sRGB
            {
                auto linear = toLinear();
                return dot(linear.rgb, getRGBLuminanceContribute(LuminanceStandard::Rec709).rgb);
            }
            else // Linear
            {
                return dot((*this).rgb, getRGBLuminanceContribute(LuminanceStandard::Rec709).rgb);
            }
        }

        /**
         * 计算与另一个颜色的对比度（Contrast Ratio）- W3C WCAG标准
         * 用于确保文本可读性
         * @param other 另一个颜色
         * @return 对比度比值 [1, 21]
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode, T>
        contrastRatio(const Color& other) const
        {
            T L1 = relativeLuminance();
            T L2 = other.relativeLuminance();
            
            // 确保L1是较亮的颜色
            if (L2 > L1) std::swap(L1, L2);
            
            // WCAG对比度公式: (L1 + 0.05) / (L2 + 0.05)
            return (L1 + T(0.05)) / (L2 + T(0.05));
        }
        

        /**
         * 将sRGB颜色转换为线性颜色空间
         * @return 线性颜色空间的Color对象
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t< rgbmode && modeIndex == 1, LinearColor<T>>
        toLinear() const
        {
            LinearColor<T> result;
            
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
            
            result.a = (*this).a; // Alpha通道保持不变
            
            return result;
        }

        /**
         * 将线性颜色空间转换为sRGB
         * @return sRGB颜色空间的Color对象
         */
        template <bool rgbmode = isRGBMode>
        typename std::enable_if_t<rgbmode && modeIndex == 0, sRGBColor<T>>
        toSRGB() const
        {
            sRGBColor<T> result;
            
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
            
            result.a = (*this).a; // Alpha通道保持不变
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
            T r = (*this).r;
            T g = (*this).g;
            T b = (*this).b;

            T max = std::max({r, g, b});
            T min = std::min({r, g, b});
            T delta = max - min;

            HSV<T> hsv;
            
            // Hue计算
            if (delta < T(1e-6))
            {
                hsv.x = T(0);
            }
            else if (max == r)
            {
                hsv.x = T(60) * (std::fmod((g - b) / delta, T(6)));
            }
            else if (max == g)
            {
                hsv.x = T(60) * ((b - r) / delta + T(2));
            }
            else
            {
                hsv.x = T(60) * ((r - g) / delta + T(4));
            }
            
            if (hsv.x < T(0))
            {
                hsv.x += T(360);
            }

            // Saturation计算 (HSV定义)
            hsv.y = (max < T(1e-6)) ? T(0) : (delta / max);

            // Value计算
            hsv.z = max;

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
            T r = (*this).r;
            T g = (*this).g;
            T b = (*this).b;

            T max = std::max({r, g, b});
            T min = std::min({r, g, b});
            T delta = max - min;

            HSL<T> hsl;
            
            // Lightness计算
            hsl.z = (max + min) / T(2);

            // Saturation计算
            if (delta < T(1e-6))
            {
                hsl.x = T(0);
                hsl.y = T(0);
            }
            else
            {
                // Hue计算（与HSV相同）
                if (max == r)
                {
                    hsl.x = T(60) * (std::fmod((g - b) / delta, T(6)));
                }
                else if (max == g)
                {
                    hsl.x = T(60) * ((b - r) / delta + T(2));
                }
                else
                {
                    hsl.x = T(60) * ((r - g) / delta + T(4));
                }
                
                if (hsl.x < T(0))
                {
                    hsl.x += T(360);
                }

                // Saturation (HSL定义)
                hsl.y = delta / (T(1) - std::abs(T(2) * hsl.z - T(1)));
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
            T r = (*this).r;
            T g = (*this).g;
            T b = (*this).b;

            HSI<T> hsi;
            
            // Intensity计算
            hsi.z = (r + g + b) / T(3);

            // Saturation计算
            T min_rgb = Min(r,g,b);
            
            if (hsi.z < T(1e-6))
            {
                hsi.y = T(0);
                hsi.x = T(0);
            }
            else
            {
                hsi.y = T(1) - min_rgb / hsi.z;
                
                // Hue计算
                if (hsi.y < T(1e-6))
                {
                    hsi.x = T(0);
                }
                else
                {
                    T numerator = T(0.5) * ((r - g) + (r - b));
                    T denominator = std::sqrt((r - g) * (r - g) + (r - b) * (g - b));
                    
                    if (denominator < T(1e-6))
                    {
                        hsi.x = T(0);
                    }
                    else
                    {
                        T theta = std::acos(numerator / denominator);
                        hsi.x = (b <= g) ? theta * T(180) / T(PI) : T(360) - theta * T(180) / T(PI);
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
                result = Vector<T, 4>{T(0.2627), T(0.6780), T(0.0593), T(1)};
                break;
            case LuminanceStandard::Rec709:
                 result = Vector<T, 4>{T(0.2126), T(0.7152), T(0.0722), T(1)};
                break;
            case LuminanceStandard::Rec601:
            default:
                result = Vector<T, 4>{T(0.299), T(0.587), T(0.114), T(1)};
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
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
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
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
        fromHSV(const Vector<T, 3>& hsv, T a = T(1))
        {
            return fromHSV(hsv.x, hsv.y, hsv.z, a);
        }

        /**
         * 从HSV Vector创建RGB颜色
         * @param hsv HSV
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
        fromHSV(const HSV<T>& hsv)
        {
            return fromHSV(hsv.x, hsv.y, hsv.z, hsv.w);
        }

        /**
         * 从HSL颜色空间创建RGB颜色
         * @param h 色相 [0, 360]
         * @param s 饱和度 [0, 1]
         * @param l 亮度 [0, 1]
         * @param a Alpha通道值 [0, 1]
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
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
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
        fromHSL(const Vector<T, 3>& hsl, T a = T(1))
        {
            return fromHSL(hsl[0], hsl[1], hsl[2], a);
        }

        /**
         * 从HSL Vector创建RGB颜色
         * @param hsl HSL
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
        fromHSL(const HSL<T>& hsl)
        {
            return fromHSL(hsl[0], hsl[1], hsl[2], hsl[3]);
        }

        /**
         * 从HSI颜色空间创建RGB颜色
         * @param h 色相 [0, 360]
         * @param s 饱和度 [0, 1]
         * @param i 强度 [0, 1]
         * @param a Alpha通道值 [0, 1]
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
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
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
        fromHSI(const Vector<T, 3>& hsi, T a = T(1))
        {
            return fromHSI(hsi[0], hsi[1], hsi[2], a);
        }

        /**
         * 从HSI Vector创建RGB颜色
         * @param hsi HSI
         * @return RGB颜色对象
         */
        template <bool rgbmode = isRGBMode>
        static typename std::enable_if_t<rgbmode, Color<T, isRGBMode, modeIndex, std::enable_if_t<std::is_arithmetic_v<T>>>>
        fromHSI(const HSI<T> hsi)
        {
            return fromHSI(hsi[0], hsi[1], hsi[2], hsi[3]);
        }

        virtual std::ostream& print(std::ostream& out) const override
        {
            std::string mode;

            if (isRGBMode)
            {
                if (modeIndex)
                {
                    mode = "sRGB";
                }
                else
                {
                    mode = "Linear";
                }
            }
            else
            {
                switch (modeIndex)
                {
                case (0):
                    mode = "HSV";
                    break;
                case (1):
                    mode = "HSL";
                    break;
                case (2):
                    mode = "HSI";
                    break;
                }
            }
			
            out << "(";
            for (size_t i = 0; i < 4; ++i) 
            {
                out << (*this).data[i];
                if (i < 4 - 1) {
                    out << ", ";
                }
            }
            out << ") (Color mode: " << mode << ")";
            return out;
        }
    };
    
}

