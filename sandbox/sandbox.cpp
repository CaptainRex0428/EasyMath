#include <iostream>
#include <ostream>

#ifdef ENGINE_API

#else
#include "EasyMath.h"

#define PRINT(var) std::cout << #var << " = " << std::endl << var << std::endl;

int main(int argc, char * argv[])
{
	using namespace EasyMath;  // 修复：使用完整命名空间名

	std::cout << "=== 测试投影矩阵函数 ===" << std::endl << std::endl;

	// ========================================
	// 测试 1: MTXLookAt - 视图矩阵
	// ========================================
	std::cout << "--- 测试 MTXLookAt (视图矩阵) ---" << std::endl;

	// 相机从侧面看原点
	Vector3 cameraPos{5.0f, 0.0f, 0.0f};      // 相机在 X 轴上
	Vector3 target{0.0f, 0.0f, 0.0f};         // 看向原点
	Vector3 up{0.0f, 1.0f, 0.0f};             // 上方向

	auto viewMatrix = MTXLookAt(cameraPos, target, up);
	PRINT(viewMatrix);

	// 验证：相机自己变换后应该在原点
	Vector<float, 4> camera_world{5.0f, 0.0f, 0.0f, 1.0f};
	auto camera_view = viewMatrix * camera_world;
	std::cout << "相机在世界空间: " << camera_world << std::endl;
	std::cout << "相机在相机空间: " << camera_view << " (应该接近原点)" << std::endl << std::endl;

	// 验证：原点在相机空间的位置
	Vector<float, 4> origin_world{0.0f, 0.0f, 0.0f, 1.0f};
	auto origin_view = viewMatrix * origin_world;
	std::cout << "原点在世界空间: " << origin_world << std::endl;
	std::cout << "原点在相机空间: " << origin_view << " (应该在相机前方)" << std::endl << std::endl;

	// ========================================
	// 测试 2: MTXPerspective - 透视投影矩阵
	// ========================================
	std::cout << "--- 测试 MTXPerspective (透视投影) ---" << std::endl;

	float fov = 3.14159f / 3.0f;     // 60 度 (PI/3)
	float aspect = 16.0f / 9.0f;     // 宽高比
	float nearPlane = 0.1f;
	float farPlane = 100.0f;

	auto perspMatrix = MTXPerspective<float>(fov, aspect, nearPlane, farPlane);  // 修复：显式指定模板参数
	PRINT(perspMatrix);

	// 测试透视效果：远处的点
	Vector<float, 4> nearPoint{1.0f, 1.0f, -1.0f, 1.0f};   // 近处
	Vector<float, 4> farPoint{1.0f, 1.0f, -50.0f, 1.0f};   // 远处

	auto nearClip = perspMatrix * nearPoint;
	auto farClip = perspMatrix * farPoint;

	std::cout << "近处点投影前: " << nearPoint << std::endl;
	std::cout << "近处点投影后: " << nearClip << std::endl;
	std::cout << "透视除法后: (" << nearClip[0]/nearClip[3] << ", "
	          << nearClip[1]/nearClip[3] << ", " << nearClip[2]/nearClip[3] << ")" << std::endl << std::endl;

	std::cout << "远处点投影前: " << farPoint << std::endl;
	std::cout << "远处点投影后: " << farClip << std::endl;
	std::cout << "透视除法后: (" << farClip[0]/farClip[3] << ", "
	          << farClip[1]/farClip[3] << ", " << farClip[2]/farClip[3] << ")" << std::endl;
	std::cout << "注意：远处的点除以 w 后坐标更小（透视效果）" << std::endl << std::endl;

	// ========================================
	// 测试 3: MTXOrtho - 正交投影矩阵
	// ========================================
	std::cout << "--- 测试 MTXOrtho (正交投影) ---" << std::endl;

	auto orthoMatrix = MTXOrtho<float>(-10.0f, 10.0f, -10.0f, 10.0f, 0.1f, 100.0f);
	PRINT(orthoMatrix);

	// 测试正交投影：远近点大小相同
	Vector<float, 4> nearPointOrtho{5.0f, 5.0f, -1.0f, 1.0f};
	Vector<float, 4> farPointOrtho{5.0f, 5.0f, -50.0f, 1.0f};

	auto nearOrtho = orthoMatrix * nearPointOrtho;
	auto farOrtho = orthoMatrix * farPointOrtho;

	std::cout << "近处点: " << nearPointOrtho << " -> " << nearOrtho << std::endl;
	std::cout << "远处点: " << farPointOrtho << " -> " << farOrtho << std::endl;
	std::cout << "注意：w 分量都是 1（无透视除法），xy 坐标相同" << std::endl << std::endl;

	// ========================================
	// 测试 4: MVP 完整变换管线
	// ========================================
	std::cout << "--- 测试完整 MVP 管线 ---" << std::endl;

	// Model: 物体变换（平移 + 旋转）
	float piOver4 = 3.14159f / 4.0f;
	auto modelMatrix = MTXTranslation(2.0f, 0.0f, -5.0f) * MTXRotationY(piOver4);

	// View: 相机在 (0, 2, 5)，看向原点
	auto view = MTXLookAt(Vector3{0.0f, 2.0f, 5.0f}, Vector3{0.0f, 0.0f, 0.0f}, Vector3{0.0f, 1.0f, 0.0f});

	// Projection: 透视投影
	float piOver3 = 3.14159f / 3.0f;
	auto proj = MTXPerspective<float>(piOver3, 16.0f / 9.0f, 0.1f, 100.0f);

	// MVP 组合
	auto MVP = proj * view * modelMatrix;

	std::cout << "Model Matrix:" << std::endl << modelMatrix << std::endl;
	std::cout << "View Matrix:" << std::endl << view << std::endl;
	std::cout << "Projection Matrix:" << std::endl << proj << std::endl;

	// 变换一个顶点
	Vector<float, 4> vertex{1.0f, 0.0f, 0.0f, 1.0f};
	auto clipSpace = MVP * vertex;

	std::cout << "顶点（物体空间）: " << vertex << std::endl;
	std::cout << "顶点（裁剪空间）: " << clipSpace << std::endl;
	std::cout << "顶点（NDC 空间）: (" << clipSpace[0]/clipSpace[3] << ", "
	          << clipSpace[1]/clipSpace[3] << ", " << clipSpace[2]/clipSpace[3] << ")" << std::endl;

	std::cout << std::endl << "=== 所有测试完成 ===" << std::endl;

	// ========================================
	// 测试 5: MTXReflection - 反射矩阵
	// ========================================
	std::cout << std::endl << "--- 测试 MTXReflection (反射矩阵) ---" << std::endl;

	// 测试 5.1: 基础反射（枚举形式）
	std::cout << ">>> 基础反射矩阵（枚举形式）：" << std::endl;

	Matrix<float, 4, 4> reflectYZ = MTXReflection<float>(ReflectionPlane::YZ);      // 翻转 X
	Matrix<float, 4, 4> reflectXZ = MTXReflection<float>(ReflectionPlane::XZ);      // 翻转 Y
	Matrix<float, 4, 4> reflectXY = MTXReflection<float>(ReflectionPlane::XY);      // 翻转 Z
	Matrix<float, 4, 4> reflectOrigin = MTXReflection<float>(ReflectionPlane::Origin);  // 翻转全部

	PRINT(reflectYZ);
	PRINT(reflectXZ);
	PRINT(reflectXY);
	PRINT(reflectOrigin);

	// 验证对角线值
	std::cout << "验证对角线值：" << std::endl;
	std::cout << "YZ 平面反射（翻转 X）：对角线 = (" << reflectYZ(0,0) << ", " << reflectYZ(1,1) << ", " << reflectYZ(2,2) << ", " << reflectYZ(3,3) << ")" << std::endl;
	std::cout << "XZ 平面反射（翻转 Y）：对角线 = (" << reflectXZ(0,0) << ", " << reflectXZ(1,1) << ", " << reflectXZ(2,2) << ", " << reflectXZ(3,3) << ")" << std::endl;
	std::cout << "XY 平面反射（翻转 Z）：对角线 = (" << reflectXY(0,0) << ", " << reflectXY(1,1) << ", " << reflectXY(2,2) << ", " << reflectXY(3,3) << ")" << std::endl;
	std::cout << "原点反射（翻转全部）：对角线 = (" << reflectOrigin(0,0) << ", " << reflectOrigin(1,1) << ", " << reflectOrigin(2,2) << ", " << reflectOrigin(3,3) << ")" << std::endl;
	std::cout << std::endl;

	// 测试 5.2: 任意平面反射
	std::cout << ">>> 任意平面反射（Householder 反射）：" << std::endl;

	Vector3 planeNormal{0.0f, 1.0f, 0.0f};      // 法向量：指向 Y 轴正方向
	Vector3 pointOnPlane{0.0f, 5.0f, 0.0f};    // 平面经过点 (0, 5, 0)，即 Y = 5 平面

	Matrix<float, 4, 4> reflectPlane = MTXReflection<float>(planeNormal, pointOnPlane);
	PRINT(reflectPlane);

	// 测试 5.3: 验证反射矩阵的性质
	std::cout << ">>> 验证反射矩阵的性质：" << std::endl;

	// 性质 1: det = -1（手性反转）
	float det = reflectPlane.determinant();
	std::cout << "行列式 det = " << det << "（应该是 -1.0）" << std::endl;

	// 性质 2: 3×3 旋转部分正交性（线性部分为 Householder 矩阵，必正交）
	// 注意：完整 4×4 反射矩阵包含平移分量，整体不是正交的——正交性只对 3×3 上三角成立
	std::cout << "3×3 旋转部分（Householder）正交性检查：" << std::endl;
	Matrix<float, 3, 3> upper3x3;
	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			upper3x3(i, j) = reflectPlane(i, j);
		}
	}
	auto orthoCheck = upper3x3 * upper3x3.transpose();
	std::cout << orthoCheck << std::endl;
	std::cout << "（3×3 上三角块应接近单位矩阵，验证线性部分的正交性）" << std::endl;

	// 性质 3: 自逆性（M == M⁻¹）
	Matrix<float, 4, 4> inverseMat = reflectPlane.inverse();
	std::cout << "自逆性检查 M == M⁻¹：" << std::endl;
	std::cout << "原矩阵：" << std::endl << reflectPlane << std::endl;
	std::cout << "逆矩阵：" << std::endl << inverseMat << std::endl;
	std::cout << "两者相等：" << (reflectPlane == inverseMat ? "是" : "否") << "（应该是 是）" << std::endl;

	// 测试 5.4: 不动点测试（平面上的点反射后不变）
	std::cout << ">>> 不动点测试（平面上的点反射后不变）：" << std::endl;

	Vector<float, 4> pointOnPlane4{0.0f, 5.0f, 0.0f, 1.0f};  // 平面上的点
	auto reflectedPoint = reflectPlane * pointOnPlane4;

	std::cout << "原点（在平面上）：" << pointOnPlane4 << std::endl;
	std::cout << "反射后：" << reflectedPoint << std::endl;
	std::cout << "验证：Y 坐标不变 = " << reflectedPoint[1] << "（应该是 5.0）" << std::endl;

	// 测试 5.5: 镜像点测试
	std::cout << ">>> 镜像点测试：" << std::endl;

	Vector<float, 4> testPoint{2.0f, 8.0f, 4.0f, 1.0f};    // 测试点 (2, 8, 4)
	auto reflectedTestPoint = reflectPlane * testPoint;

	std::cout << "原点：" << testPoint << std::endl;
	std::cout << "反射后：" << reflectedTestPoint << std::endl;
	std::cout << "验证：X 坐标不变 = " << reflectedTestPoint[0] << "（应该是 2.0）" << std::endl;
	std::cout << "验证：Y 坐标镜像 = " << reflectedTestPoint[1] << "（应该是 2.0 = 5 - (8-5)）" << std::endl;
	std::cout << "验证：Z 坐标不变 = " << reflectedTestPoint[2] << "（应该是 4.0）" << std::endl;

	// 测试 5.6: 一致性测试（枚举形式 vs 通用形式）
	std::cout << ">>> 一致性测试（枚举形式 vs 通用形式）：" << std::endl;

	// YZ 平面反射（枚举）vs 法向量 (-1, 0, 0) 过原点的平面（通用）
	Matrix<float, 4, 4> reflectYZ_enum = MTXReflection<float>(ReflectionPlane::YZ);
	Vector3 yzNormal{-1.0f, 0.0f, 0.0f};
	Vector3 yzPoint{0.0f, 0.0f, 0.0f};
	Matrix<float, 4, 4> reflectYZ_generic = MTXReflection<float>(yzNormal, yzPoint);

	std::cout << "枚举形式 MTXReflection(YZ)：" << std::endl << reflectYZ_enum << std::endl;
	std::cout << "通用形式 MTXReflection({-1,0,0}, {0,0,0})：" << std::endl << reflectYZ_generic << std::endl;
	std::cout << "两者相等：" << (reflectYZ_enum == reflectYZ_generic ? "是" : "否") << "（应该是 是）" << std::endl;

	return 0;
}

#endif
