# 使用指南

## 只计算一条谱的 Q：下载一个文件即可

- MATLAB：下载根目录 [qcalc_gui.m](../qcalc_gui.m)，在 MATLAB 中打开后运行。建议 MATLAB R2020b+，无需额外工具箱。
- Python：下载根目录 [qcalc_gui.py](../qcalc_gui.py)，运行 `python qcalc_gui.py`。需要 NumPy、SciPy、Matplotlib、Tkinter；Windows 常规 Python/Anaconda 一般自带 Tkinter。Linux 图形界面还需要可用的显示环境。

两个文件分别独立运行，**不需要同时下载，也不需要其他仓库源码文件**。Python 首次安装依赖可运行：

```bash
python -m pip install numpy scipy matplotlib
python qcalc_gui.py
```

## 操作

1. 打开你的单组 TXT/CSV 数据。`p260r110.txt` 的默认设置为幅度、THz、频率列 1、数据列 2。
2. 如果文件实际为功率、dB 或复数，选择对应类型和列号后重新读取。内部按功率谱拟合：幅度平方、dB 转 `10^(dB/10)`。
3. 输入频率窗口上下限，单位 THz。Python 也可以在图上拖动选区。
4. 多峰时点选目标峰并设定“选中峰邻域”，或手动缩小窗口，避免把多个模式拟合成单个共振。
5. 选择计算方式；自定义拟合可选择 Fano、Lorentzian、Gaussian、Voigt、pseudo-Voigt 和常数/线性背景。
6. 计算后对照拟合曲线和残差，滚动查看定义与可靠性提示。稍微调整窗口，检查 Q 是否稳定。
7. 导出 JSON、CSV、PNG；输入文件由你自行提供，不随仓库发布。

Fano 的极点 Q=`f0/Gamma`；半突出度、拟合半对比度、绝对半功率 Q 是不同指标，不能混用。Python 的“筛选同算法”按原局部窗口规则计算；MATLAB 和 Python 自定义拟合使用整个指定区间，不保证不同窗口/优化器给出相同数值。导出点密不代表求解器的真实采样足够密。

## 批处理、函数接口与无图形服务器

原有接口仍保留，下载整个仓库后使用：

```bash
python fit_q.py /path/to/S21.txt --kind amplitude --model fano --feature peak --window 0.468 0.482 --output results/sample --no-plot
python python/resonance_qfit/batch_qfit.py /path/to/spectra --pattern "*.txt" --kind amplitude --model fano --window 0.468 0.482 --output-dir results/batch --no-plots
```

窗口只是示例，请改成目标共振区间。命令行版本不需要 Tkinter 图形显示；Python 依赖见 `python/resonance_qfit/requirements.txt`。

MATLAB 批处理和直接调用见 [MATLAB 工具说明](../matlab/resonance_qfit/README.md)。两份偏振数据自动比较保留在根目录 `EXAMPLE_P260R110_TWO_POL.m`，需要完整仓库及你自己的两份数据。

## 其他功能

根目录 `START_HERE.m` 仍提供任务菜单；Q 的交互入口统一指向新版单文件界面。

| 功能 | 位置 |
|---|---|
| THz 时域模拟信号 | `matlab/tds/data_sin_shengcheng.m` |
| THz-TDS 数据处理 | `matlab/tds/shujuchuli.m` |
| CST 通用参数采样与扫描模板 | `matlab/cst/` |
| 三开口环槽采样与扫描 | [matlab/cst/sankaiKou/](../matlab/cst/sankaiKou/README.md) |
| 模型示意图绘制 | `visualization/` |

通用 CST 模板的参数与结果树路径需要按实际工程修改；三开口环槽请优先参考专用实例。TDS 处理同样需设置实际列号、单位、FFT 区间及是否加窗。

删除的旧 GUI 向导、重复单谱入口和固定路径单谱示例均已被 `qcalc_gui.m` 替代，历史版本可在 Git 提交记录中找回。数值测试、批处理、函数库及不同任务工具保留。

Q 定义、公式和限制见 [详细文档](resonance_qfit.md)。
