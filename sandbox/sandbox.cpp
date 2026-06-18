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

	return 0;
}

#endif
