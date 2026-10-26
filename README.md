# THz-Metasurface-AI

太赫兹超表面研究代码仓库。目前收录 THz-TDS 模拟信号生成和时域数据处理的 MATLAB 脚本，以及基于 Blender 的结构示意图建模与渲染脚本。Python 正向/逆向模型尚未上传。

## 目录

- `matlab/tds/data_sin_shengcheng.m`：生成模拟参考/样品时域信号，绘制时域与频域图，并在 MATLAB 当前工作目录写出 `THz_simulated.xlsx`。
- `matlab/tds/shujuchuli.m`：读取时域数据，检查采样与单位，截取信号、可选加高斯窗、计算 FFT 幅度谱和幅度透射谱。
- `matlab/cst/`：通用参数生成与 MATLAB–CST 联动扫描脚本；配置步骤和 AuSRR 随机性评估见 [`matlab/cst/README.md`](matlab/cst/README.md)。
- `python/resonance_qfit/`：Python 共振线形拟合、Q 值计算和批处理，支持 Lorentzian、Fano、Gaussian、Voigt、pseudo-Voigt。
- `matlab/resonance_qfit/`：基础 MATLAB 实现的同类 Q 值工具，无需额外拟合工具箱；运行 `run_qfit` 可选择 CST 导出文件。
- `visualization/render_metasurface.py`：使用 Blender 生成 4×4 阵列示意模型、渲染图和可编辑模型文件；运行方法见 `visualization/README.md`。

## 使用 MATLAB 示例

1. 在 MATLAB 中运行 `matlab/tds/data_sin_shengcheng.m`。生成的 Excel 文件是演示数据，不纳入版本控制。
2. 若要用 `shujuchuli.m` 处理该三列文件，请先在脚本顶部设置 `timeCol = 1`、`refCol = 2`、`samCol = 3`、`fftStart = 0`、`fftEnd = 20e-12`、`useGaussianWindow = false`，然后运行脚本并选择刚生成的文件。
3. 处理实测数据时，按实际文件调整列号、时间/信号单位和 FFT 区间；如开启高斯窗，也要检查 `gaussianCenter` 与数据时间范围是否相符。脚本默认值针对原作者的实测数据，不能直接用于上述三列演示文件。

脚本使用 MATLAB 的表格读写、绘图和 FFT 功能；尚未针对具体 MATLAB 版本进行运行验证。

## 数据与凭据

共振 Q 工具的输入格式、公式和用法见 [共振线形与 Q 值说明](docs/resonance_qfit.md)。特别注意输入幅值与功率的区别，以及 Fano 的极点线宽 Q 与可见峰/谷宽度 Q 的区别。

`.gitignore` 排除常见原始/生成数据、模型权重、CST 工程文件、凭据和本地环境。添加新文件前仍需逐项检查内容；忽略规则不会阻止显式强制提交。
