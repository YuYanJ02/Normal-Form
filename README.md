# catalogue

地月圆型限制性三体问题（CR3BP）平动点附近哈密顿量的 **Lie 变换正规形** 与 **作用-角变量** 计算代码，包含 MATLAB 符号/数值流程与 C++ 高性能实现，可用于周期/拟周期轨道的动力学分析与 catalogue 构建。

## 功能概览

- **哈密顿量展开**：在平动点中心坐标系下，对 CRTBP 哈密顿量进行高阶多项式展开。
- **对角化与复数化**：对线性部分对角化，并在复数坐标下表示高阶项。
- **Lie 正规形**：基于 Dragt–Finn Lie 级数计算生成函数，提取共振项得到正规形哈密顿量。
- **作用-角变量变换**：从复数正规形坐标构造作用-角变量，并得到 \(H(I)\) 形式。
- **坐标变换链**：共转坐标 → 正则坐标 → 对角化坐标 → 复数坐标 → 正规形坐标 → 作用-角变量。
- **C++ 实现**：多项式表示、Poisson 括号、Lie 变换与正规形的 C++ 实现，可单独编译运行或与 MATLAB 交互。

## 目录结构

```text
catalogue/
├── main.m                          # 主流程：哈密顿量 → 复数 → 正规形 → 作用-角
├── main_Orbit_.m                   # 轨道积分，生成 qp 采样点
├── compute_libration_point_params.m
├── expand_hamiltonian.m
├── hamiltonian_after_diagonalization.m
├── complexify_hamiltonian.m
├── lie_transformation.m            # MATLAB 版 Lie 级数正规形
├── transform_to_action_angle.m
├── build_coordinate_transforms.m
├── crtbp_state_to_complex.m
├── fun_crtbp_.m                    # CRTBP 动力学方程
├── export_H_complex.m
├── reduce.m, plot_action_angle.m
├── action_angle.txt, H_complex_terms.txt, qp_points.txt
├── CMakeLists.txt
├── .vscode/                        # VS Code 调试配置
└── Normal_Form/                    # C++/Visual Studio 正规形工程
    ├── normal_form/                # 核心 C++ 源码
    │   ├── main.cpp
    │   ├── normal_form.cpp
    │   └── normal_form.hpp
    └── normal_form.sln             # VS 工程文件
```

## 使用说明

### 1. MATLAB 正规形流程

1. 在 MATLAB 中将 `catalogue` 目录加入路径。
2. 打开并运行 `main.m`：
   - 计算平动点参数、哈密顿量展开与对角化；
   - 进行复数化与（可选的）Lie 正规形变换；
   - 生成 `H_complex_terms.txt`、`qp_points.txt` 等中间结果。
3. 如需生成轨道采样点，可运行 `main_Orbit_.m`，再用 `crtbp_state_to_complex.m` 将轨道状态转换为复数坐标。

### 2. C++ 正规形模块

1. 安装支持 C++17 的编译器（MSVC / GCC / Clang）。
2. 使用 CMake：在命令行中进入 `catalogue`，执行：
   ```bash
   cmake -S . -B build
   cmake --build build --config Release
   ```
   或直接用 `Normal_Form/normal_form.sln` 在 Visual Studio 中打开并编译。
3. 运行生成的可执行程序，可读入 MATLAB 导出的 `H_complex_terms.txt` / `qp_points.txt`，完成正规形与作用-角变量的数值计算。

## 依赖环境

- MATLAB（建议 R2019a 及以上；如使用符号推导需要 Symbolic Math Toolbox）。
- C++17 编译器与 CMake（可选，仅用于 C++ 正规形部分）。

## 适用场景

适用于对 **地月 CR3BP 平动点邻域** 的精细动力学建模：周期轨道 catalogue 构建、长期稳定性分析、共振结构研究等，也可作为正规形与 Lie 变换在航天动力学中的教学示例。
