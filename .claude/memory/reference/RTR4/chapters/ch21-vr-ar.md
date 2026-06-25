# Ch21 Virtual and Augmented Reality VR/AR

**核心问题**：VR/AR 渲染有什么特殊挑战？

---

## 21.1 设备和系统概述

**VR 类型**：
- **PC VR**：高保真，需要高性能 PC
- **Standalone VR**：Quest 等，自带 GPU
- **手机 VR**：Cardboard 等，已过时

**AR 类型**：
- 手机/平板
- AR 眼镜（HoloLens、Magic Leap）
- 透明度 + 真实世界合成

---

## 21.2 物理元素

### 21.2.1 延迟（Latency）

**核心问题**：运动到光子（Motion-to-Photon）的总延迟。

**目标**：< 20ms（最好 < 11ms）。

**延迟源**：
- 传感器采样
- 姿态预测
- CPU 渲染
- GPU 渲染
- 扫描输出
- 显示器响应

**缓解**：
- **异步时间扭曲**（Asynchronous Timewarp, ATW）
- **异步空间扭曲**（ASW）
- **运动预测**

### 21.2.2 光学（Optics）

**焦距-调节冲突**（Vergence-Accommodation Conflict, VAC）：
- 双眼聚焦于虚拟物体（近）
- 晶状体却适应远（无穷远）
- 导致眼睛疲劳

**缓解**：
- 变焦显示
- 光场显示
- 多焦平面

**视场**（FOV）：约 100-120°（小于人眼 ~200°，但高于普通显示器）

### 21.2.3 立体视觉（Stereo Vision）

**核心**：左眼右眼各看一张图 → 立体感。

**视差**（Parallax）：左右眼之间的位置差。

**会聚距离**（Vergence Distance）：两眼视线的交点。

---

## 21.3 API 和硬件

### 21.3.1 立体渲染

**传统**：
- 渲染两次（左眼、右眼）
- 各有自己的 view/proj 矩阵
- 时间扭曲校正

**Single-Pass Stereo**：
- 一次提交渲染两个视点
- 用 SV_ViewportArrayIndex 或实例化

### 21.3.2 注视点渲染（Foveated Rendering）

**核心**：根据眼睛注视点分配渲染分辨率。

**注视区域**：高分辨率
**周边**：低分辨率

**实现**：
- 眼动追踪
- 多分辨率渲染（tiled 或 shifted）
- 节省 30-50% GPU

---

## 21.4 渲染技术

### 21.4.1 抖动（Dithering）

**核心**：用抖动掩饰低分辨率（如 6DoF 跟踪不足）。

- 8 色抖动
- Bayer 矩阵
- 配合时间累积

### 21.4.2 计时（Timing）

**预测**：
- 渲染提前开始（用预测姿态）
- 渲染时姿态稍变 → 用 ATW 校正

**时间扭曲**（Time Warp）：
- 渲染完成后，扫描输出前
- 重新投影到最新姿态
- 仅做姿态修正（不含场景变化）

**空间扭曲**（Space Warp）：
- 推断运动场景（基于深度+姿态）
- 补偿物体运动

---

## 与 EasyMath 关联

**低**——这是 VR 渲染层。

**但是**：
- 立体渲染需要左右眼 View/Proj 矩阵
- EasyMath 可加左右眼矩阵生成

## 关键术语

| 术语 | 释义 |
|------|------|
| VR | Virtual Reality |
| AR | Augmented Reality |
| 延迟 | 运动到光子时间 |
| ATW | 异步时间扭曲 |
| ASW | 异步空间扭曲 |
| VAC | 焦距-调节冲突 |
| 立体视觉 | 左右眼不同图 |
| 视差 | 双眼位置差 |
| 注视点渲染 | Foveated Rendering |
| 时间扭曲 | Time Warp |
| 空间扭曲 | Space Warp |
