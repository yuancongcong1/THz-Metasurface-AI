# THz-Metasurface-AI

太赫兹超表面研究代码仓库。目前收录 THz-TDS 模拟信号生成和时域数据处理的 MATLAB 脚本。Python 正向/逆向模型与 CST 自动化脚本尚未上传。

## 目录

- `matlab/tds/data_sin_shengcheng.m`：生成模拟参考/样品时域信号，绘制时域与频域图，并在 MATLAB 当前工作目录写出 `THz_simulated.xlsx`。
- `matlab/tds/shujuchuli.m`：读取时域数据，检查采样与单位，截取信号、可选加高斯窗、计算 FFT 幅度谱和幅度透射谱。

## 使用 MATLAB 示例

1. 在 MATLAB 中运行 `matlab/tds/data_sin_shengcheng.m`。生成的 Excel 文件是演示数据，不纳入版本控制。
2. 若要用 `shujuchuli.m` 处理该三列文件，请先在脚本顶部设置 `timeCol = 1`、`refCol = 2`、`samCol = 3`、`fftStart = 0`、`fftEnd = 20e-12`、`useGaussianWindow = false`，然后运行脚本并选择刚生成的文件。
3. 处理实测数据时，按实际文件调整列号、时间/信号单位和 FFT 区间；如开启高斯窗，也要检查 `gaussianCenter` 与数据时间范围是否相符。脚本默认值针对原作者的实测数据，不能直接用于上述三列演示文件。

脚本使用 MATLAB 的表格读写、绘图和 FFT 功能；尚未针对具体 MATLAB 版本进行运行验证。

## 数据与凭据

`.gitignore` 排除常见原始/生成数据、模型权重、CST 工程文件、凭据和本地环境。添加新文件前仍需逐项检查内容；忽略规则不会阻止显式强制提交。
