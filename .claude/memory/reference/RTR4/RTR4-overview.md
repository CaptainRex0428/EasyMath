# RTR4 (Real-Time Rendering, 4th Edition) 总览

> **源文件**：`reference/Real-Time Rendering, Fourth Edition (Tomas Akenine-Moller, Eric Haines, Naty Hoffman etc.).epub`
> **作者**：Tomas Akenine-Möller, Eric Haines, Naty Hoffman, Angelo Pesce, Michal Iwanicki, Sébastien Hillaire
> **中文译本**：实时渲染第四版
> **生成方式**：pandoc 转换（中文版）
> **状态**：本 memory 已抽取概要，**未来提取相关知识点时不再读取源文件**
> **章节数**：26 章正文 + 序言/标题页

---

## 用途说明（与 EasyMath 项目的关联）

RTR4 是实时渲染领域的**核心参考书**。当本项目（EasyMath）相关的任务涉及以下领域时，**先查阅 `.claude/memory/reference/RTR4/`**，而不是回读源文件：

| 任务领域 | 主要参考章节 |
|---------|-------------|
| 矩阵/变换/四元数 实现与 API 设计 | Ch4 Transform, Ch2 管线 |
| 投影矩阵 / 视图矩阵推导 | Ch4 §4.7 |
| 颜色空间 / 光度学 / 亮度 | Ch8 Light and Color |
| 物理着色（BRDF、菲涅尔、微表面） | Ch9 Physically Based Shading |
| 光照模型 | Ch10 局部 / Ch11 全局 |
| 阴影算法 | Ch7 Shadows |
| 纹理过滤（mipmap、各向异性） | Ch6 Texturing |
| 加速结构 / 空间分割 | Ch19 Acceleration, Ch22 Intersection |
| LOD / 剔除 | Ch19.9 / Ch19.2-7 |
| 碰撞检测 | Ch25 |
| 光线追踪 | Ch26 |
| 性能分析 | Ch18 |
| GPU 架构 | Ch23 |

---

## 全书章节快速索引（26 章）

### 第一部分：基础（Foundations）

| Ch | 英文标题 | 中文标题 | 与 EasyMath 关联度 |
|----|---------|---------|------------------|
| 1 | Introduction | 简介 | 低（铺垫） |
| 2 | The Graphics Rendering Pipeline | 图形渲染管线 | 中（架构总览） |
| 3 | The Graphics Processing Unit | 图形处理单元 | 低（硬件） |
| 4 | **Transform** | **变换** | **★ 极高**（与 Matrix.h / Quaternion.h 强相关） |
| 5 | Shading Basics | 着色基础 | 中 |
| 6 | Texturing | 纹理 | 中 |
| 7 | Shadows | 阴影 | 中（不含代码） |

### 第二部分：光与颜色（Light and Color）

| Ch | 英文标题 | 中文标题 | 与 EasyMath 关联度 |
|----|---------|---------|------------------|
| 8 | **Light and Color** | **光与颜色** | **★ 高**（与 Color.h 相关） |
| 9 | **Physically Based Shading** | **基于物理的着色** | **★ 高**（BRDF 公式参考） |
| 10 | Local Illumination | 局部光照 | 中 |
| 11 | Global Illumination | 全局光照 | 中 |
| 12 | Image-Space Effects | 图像空间特效 | 低 |

### 第三部分：几何与处理（Geometry and Processing）

| Ch | 英文标题 | 中文标题 | 与 EasyMath 关联度 |
|----|---------|---------|------------------|
| 13 | Beyond Polygons | 超越多边形 | 低（粒子/体素） |
| 14 | Volumetric and Translucency | 体积与半透明 | 低 |
| 15 | Non-Photorealistic Rendering | 非真实感渲染 | 低 |
| 16 | Polygonal Techniques | 多边形技术 | 低 |
| 17 | Curves and Curved Surfaces | 曲线和曲面 | 中（Bezier/B-Spline） |

### 第四部分：高级与优化（Advanced and Optimization）

| Ch | 英文标题 | 中文标题 | 与 EasyMath 关联度 |
|----|---------|---------|------------------|
| 18 | Pipeline Optimization | 管线优化 | 中 |
| 19 | **Acceleration Algorithms** | **加速算法** | **★ 高**（BVH/kd-tree/八叉树/LOD） |
| 20 | Efficient Shading | 高效着色 | 中 |
| 21 | Virtual and Augmented Reality | VR/AR | 低 |
| 22 | **Intersection Test Methods** | **相交测试方法** | **★ 高**（射线/包围体/平面） |
| 23 | Graphics Hardware | 图形硬件 | 低 |
| 24 | The Future | 未来 | 低（展望） |
| 25 | **Collision Detection** | **碰撞检测** | **★ 中**（GJK/SAT/BVH CD） |
| 26 | **Real-Time Ray Tracing** | **实时光线追踪** | **★ 中**（BVH 遍历/Coherence） |

---

## 章节详细索引

### Ch1 Introduction 简介
- 1.1 内容概述
- 1.2 符号和定义
  - 1.2.1 数学符号
  - 1.2.2 几何定义
  - 1.2.3 着色
- 深入阅读和资源

### Ch2 The Graphics Rendering Pipeline 图形渲染管线
- 2.1 渲染管线的架构
- 2.2 应用阶段
- 2.3 几何处理阶段
  - 2.3.1 顶点着色
  - 2.3.2 可选的顶点处理
  - 2.3.3 裁剪
  - 2.3.4 屏幕映射
- 2.4 光栅化阶段
  - 2.4.1 三角形设置
  - 2.4.2 三角形遍历
- 2.5 像素处理阶段
  - 2.5.1 像素着色
  - 2.5.2 合并
- 2.6 回顾整个管线
- 2.7 总结
- 2.8 补充阅读和资源

### Ch3 The Graphics Processing Unit 图形处理单元
- 3.1 数据并行结构
- 3.2 GPU 管线概述
- 3.3 可编程着色器阶段
- 3.4 可编程着色及其 API 的演变
- 3.5 顶点着色器
- 3.6 曲面细分阶段
- 3.7 几何着色器
  - 3.7.1 流式输出
- 3.8 像素着色器
- 3.9 合并阶段
- 3.10 计算着色器

### Ch4 Transform 变换（**核心参考**）
- 4.1 基本变换
  - 4.1.1 平移
  - 4.1.2 旋转
  - 4.1.3 缩放
  - 4.1.4 剪切
  - 4.1.5 变换的连接
  - 4.1.6 刚体变换
  - 4.1.7 法线变换
  - 4.1.8 计算逆矩阵
- 4.2 特殊的矩阵变换和操作
  - 4.2.1 欧拉变换
  - 4.2.2 从欧拉变换中提取参数
  - 4.2.3 矩阵分解
  - 4.2.4 绕任意轴旋转
- 4.3 四元数
  - 4.3.1 数学背景
  - 4.3.2 四元数变换（矩阵转换、slerp、向量旋转）
- 4.4 顶点混合
- 4.5 变形
- 4.6 几何缓存回放
- 4.7 投影
  - 4.7.1 正交投影
  - 4.7.2 透视投影

### Ch5 Shading Basics 着色基础
- 5.1 着色模型
- 5.2 光源（方向光、点光源、聚光灯）
- 5.3 实现着色模型（计算频率、材质系统）
- 5.4 锯齿和抗锯齿（采样滤波理论、屏幕空间 AA、形态学方法）
- 5.5 透明度、Alpha、合成（混合顺序、OIT、Alpha 预乘）
- 5.6 显示编码

### Ch6 Texturing 纹理
- 6.1 纹理管线（投影函数、转换函数、纹理值）
- 6.2 图像纹理（放大、缩小、Mipmap、SAT、各向异性、立方体贴图、压缩）
- 6.3 程序化纹理
- 6.4 纹理动画
- 6.5 材质映射
- 6.6 Alpha 映射
- 6.7 凹凸映射（Blinn 法线映射）
- 6.8 视差映射
- 6.9 纹理光源

### Ch7 Shadows 阴影
- 7.1 平面阴影（投影阴影、软阴影）
- 7.2 曲面上的阴影
- 7.3 阴影体算法
- 7.4 阴影贴图
- 7.5 PCF
- 7.6 PCSS
- 7.7 过滤阴影贴图
- 7.8 体积阴影
- 7.9 不规则 z-buffer
- 7.10 其他应用

### Ch8 Light and Color 光与颜色（**核心参考**）
- 8.1 光量
  - 8.1.1 辐射度量学（Radiometry）
  - 8.1.2 光度学（Photometry）
  - 8.1.3 色度学（Colorimetry）
  - 8.1.4 使用 RGB 颜色进行渲染
- 8.2 从场景到屏幕
  - 8.2.1 HDR 显示编码
  - 8.2.2 色调映射
  - 8.2.3 颜色分级

### Ch9 Physically Based Shading 基于物理的着色（**核心参考**）
- 9.1 光的物理学（粒子、介质、表面、次表面散射）
- 9.2 相机
- 9.3 The BRDF
- 9.4 光照
- 9.5 菲涅尔反射
  - 9.5.1 外反射
  - 9.5.2 典型菲涅尔值（电介质/金属/半导体/水中/参数化）
  - 9.5.3 内反射
- 9.6 微观几何
- 9.7 微表面理论
- 9.8 表面反射的 BRDF 模型
  - 9.8.1 法线分布函数（NDF，各向同性/各向异性）
  - 9.8.2 多次反弹
- 9.9 次表面散射 BRDF
- 9.10 布料 BRDF（经验/微表面/微圆柱体）
- 9.11 波动光学 BRDF（衍射/薄膜干涉）
- 9.12 分层材质
- 9.13 混合和过滤材质

### Ch10 Local Illumination 局部光照
- 10.1 面光源（光泽材质/一般形状）
- 10.2 环境光照
- 10.3 球面/半球函数
  - 10.3.1 简单表格
  - 10.3.2 球面基底（球面径向基/球面高斯/球谐）
  - 10.3.3 半球基底（AHD/HL2/H-Basis）
- 10.4 环境映射（经纬度/球面/立方体/其他）
- 10.5 基于图像的高光照明（IBL，预过滤环境贴图，BRDF 分裂积分近似）
- 10.6 irradiance 环境映射（球谐 IBL）
- 10.7 误差来源

### Ch11 Global Illumination 全局光照
- 11.1 渲染方程
- 11.2 通用 GI（辐射度、光线追踪）
- 11.3 环境光遮蔽（SSAO、HBAO 等）
- 11.4 定向遮蔽
- 11.5 漫反射 GI（SH/LPV/Voxel/Screen-space）
- 11.6 镜面 GI（环境贴图反射、SSR、平面反射）
- 11.7 统一方法

### Ch12 Image-Space Effects 图像空间特效
- 12.1 图像处理（双边滤波）
- 12.2 重投影技术
- 12.3 镜头光晕和泛光
- 12.4 景深
- 12.5 运动模糊

### Ch13 Beyond Polygons 超越多边形
- 13.1 渲染频谱
- 13.2 固定视图效果
- 13.3 天空盒
- 13.4 光场渲染
- 13.5 Sprite 和图层
- 13.6 广告牌技术（屏幕对齐/世界朝向/轴向/Impostor）
- 13.7 位移技术
- 13.8 粒子系统（着色、模拟）
- 13.9 点渲染
- 13.10 体素

### Ch14 Volumetric and Translucency 体积与半透明
- 14.1 光线散射理论（参与介质、透光率、散射、相位函数——瑞利/米氏/几何）
- 14.2 特殊体渲染（雾、体积光）
- 14.3 通用体渲染
- 14.4 天空渲染（空气透视、云）
- 14.5 半透明表面（折射、焦散）
- 14.6 次表面散射（环绕光照、法线模糊、预积分皮肤、纹理空间扩散、屏幕空间扩散）
- 14.7 毛发和皮毛
- 14.8 统一方法

### Ch15 Non-Photorealistic Rendering 非真实感渲染
- 15.1 卡通着色
- 15.2 轮廓渲染（基于法线/几何/图像处理边缘检测）
- 15.3 笔触表面风格化
- 15.4 线条
- 15.5 文本渲染

### Ch16 Polygonal Techniques 多边形技术
- 16.1 三维数据来源
- 16.2 曲面细分和三角形划分（着色问题、边界开裂、T 顶点）
- 16.3 整合（合并、定向、实体性、法线平滑）
- 16.4 三角形扇/带/网格
- 16.5 简化
- 16.6 压缩和精度

### Ch17 Curves and Curved Surfaces 曲线和曲面
- 17.1 参数化曲线
  - 17.1.1 Bezier 曲线（Bernstein/有理 Bezier）
  - 17.1.2 GPU 上的有界 Bezier
  - 17.1.3 连续性与分段 Bezier
  - 17.1.4 三次 Hermite 插值
  - 17.1.5 Kochanek-Bartels 曲线
  - 17.1.6 B-样条
- 17.2 参数化曲面（Bezier 面片/三角形、PN 三角形、Phong 细分、B-样条曲面）
- 17.3 隐式表面
- 17.4 细分曲线
- 17.5 细分表面（Loop/Catmull-Clark/分段平滑/位移细分）
- 17.6 高效曲面细分（分数细分/自适应/快速 Catmull-Clark）

### Ch18 Pipeline Optimization 管线优化
- 18.1 分析和调试工具
- 18.2 定位性能瓶颈（按阶段测试）
- 18.3 性能测量
- 18.4 优化（应用/API/几何/光栅化/像素/帧缓冲/合并阶段）
- 18.5 多处理（流水线/并行/基于任务/API 支持）

### Ch19 Acceleration Algorithms 加速算法（**核心参考**）
- 19.1 空间数据结构
  - 19.1.1 层次包围体
  - 19.1.2 BSP 树（k-d 树/多边形对齐 BSP）
  - 19.1.3 八叉树
  - 19.1.4 缓存无关/缓存感知表示
  - 19.1.5 场景图
- 19.2 剔除技术
- 19.3 背面剔除
- 19.4 视锥体剔除
- 19.5 入口剔除
- 19.6 细节剔除和小三角形剔除
- 19.7 遮挡剔除（遮挡查询/层次 Z 缓冲 HiZ）
- 19.8 剔除系统
- 19.9 LOD（离散/混合/Alpha/CLOD/地貌；选择：基于范围/投影面积）
- 19.10 渲染大型场景（虚拟纹理、流式传输、地形）

### Ch20 Efficient Shading 高效着色
- 20.1 延迟着色
- 20.2 贴花渲染
- 20.3 分块着色（Tiled Shading）
- 20.4 聚类着色（Clustered Shading）
- 20.5 延迟纹理
- 20.6 对象空间和纹理空间着色

### Ch21 Virtual and Augmented Reality VR/AR
- 21.1 设备和系统概述
- 21.2 物理元素（延迟、光学、立体视觉）
- 21.3 API 和硬件（立体渲染、注视点渲染）
- 21.4 渲染技术（抖动、计时）

### Ch22 Intersection Test Methods 相交测试方法（**核心参考**）
- 22.1 GPU 加速的拾取
- 22.2 定义和工具
- 22.3 创建包围体（AABB/k-DOP/球/凸多面体/OBB）
- 22.4 几何概率
- 22.5 经验法则
- 22.6 射线/球体相交（数学解法/优化解法）
- 22.7 射线/Box 相交（平板法/斜率法）
- 22.8 射线/三角形相交（Möller-Trumbore 等）
- 22.9 射线/多边形相交
- 22.10 平面/Box 相交（AABB/OBB）
- 22.11 三角形/三角形相交
- 22.12 三角形/Box 相交
- 22.13 BV/BV 相交（球-球、球-Box、AABB-AABB、k-DOP、OBB-OBB）
- 22.14 视锥体相交测试
- 22.15 线/线相交（2D/3D）
- 22.16 三平面相交

### Ch23 Graphics Hardware 图形硬件
- 23.1 光栅化（插值/保守光栅化）
- 23.2 大规模计算和调度
- 23.3 延迟和占用率
- 23.4 内存架构和总线
- 23.5 缓存和压缩
- 23.6 颜色缓冲（视频显示控制器、单/双/三重缓冲）
- 23.7 深度剔除、测试和缓冲
- 23.8 纹理化
- 23.9 架构
- 23.10 案例分析（ARM Mali G71 / NVIDIA Pascal / AMD GCN Vega）
- 23.11 光线追踪架构

### Ch24 The Future 未来
- 24.1 其他事项
- 24.2 你

### Ch25 Collision Detection 碰撞检测（**核心参考**）
- 25.1 宽阶段（扫描剪枝、网格、层次包围体）
- 25.2 中阶段（BVH 构建、BVH 测试、成本函数、OBB 树）
- 25.3 窄阶段（图元 vs 图元、距离查询）
- 25.4 射线碰撞检测
- 25.5 BSP 树动态 CD
- 25.6 限时碰撞检测
- 25.7 可变形模型
- 25.8 连续碰撞检测
- 25.9 碰撞响应
- 25.10 粒子（系统、物理模拟）
- 25.11 动态相交测试（球-平面、球-球、球-多边形、SAT）

### Ch26 Real-Time Ray Tracing 实时光线追踪
- 26.1 光线追踪基础
- 26.2 光线追踪着色器
- 26.3 顶层和底层加速结构
- 26.4 一致性（场景一致性/空间数据结构属性/构造方案/遍历方案）

---

## 配套文件

- `chapters/ch01-...md` ~ `chapters/ch26-...md`：每章详细概要
- `formulas.md`：跨章节关键公式集合（LaTeX + 来源章节 + C++ 实现片段）
- `principles.md`：核心原理/算法思想集合

## 提取完成度

- [x] Ch1 (TOC)
- [x] Ch2 (TOC)
- [x] Ch3 (TOC)
- [x] Ch4 (TOC, 详尽版待补)
- [x] Ch5 (TOC)
- [x] Ch6 (TOC)
- [x] Ch7 (TOC)
- [x] Ch8 (TOC)
- [x] Ch9 (TOC)
- [x] Ch10 (TOC)
- [x] Ch11 (TOC)
- [x] Ch12 (TOC)
- [x] Ch13 (TOC)
- [x] Ch14 (TOC)
- [x] Ch15 (TOC)
- [x] Ch16 (TOC)
- [x] Ch17 (TOC)
- [x] Ch18 (TOC)
- [x] Ch19 (TOC)
- [x] Ch20 (TOC)
- [x] Ch21 (TOC)
- [x] Ch22 (TOC)
- [x] Ch23 (TOC)
- [x] Ch24 (TOC)
- [x] Ch25 (TOC)
- [x] Ch26 (TOC)
- [ ] 各章详细概要（基于源文件）
- [ ] formulas.md 关键公式
- [ ] principles.md 核心原理
