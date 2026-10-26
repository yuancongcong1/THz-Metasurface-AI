# MATLAB 共振 Q 值工具

新手建议直接打开仓库根目录 `START_QFIT.m` 点“运行”，按弹窗选择文件和数据类型，在频谱上选择拟合窗口。`EXAMPLE_P260R110.m` 已设置这份示例的 0.466–0.471 THz 窗口；输入文件需自行提供。

逐步说明见 [第一次使用指南](../../docs/START_HERE.md)。

在仓库根目录的 MATLAB 命令窗口运行：

```matlab
addpath('matlab/resonance_qfit');
run_qfit
```

在 `run_qfit.m` 顶部选择幅值/功率/dB/复数和线形。支持 Lorentzian、Fano、Gaussian、Voigt、pseudo-Voigt。`fit_resonance_q.m` 可直接处理数组，`batch_qfit.m` 可批处理。

不依赖额外拟合工具箱；使用基础 MATLAB `fminsearch`，多初值优化。演示：`demo_qfit`；数值和读写测试：`test_qfit`。

Fano 的拟合极点 Q 与可见峰/谷宽度 Q 可能不同。完整用法、公式、输出定义和限制见 [说明文档](../../docs/resonance_qfit.md)。
