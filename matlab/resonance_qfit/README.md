# MATLAB 共振 Q 值工具

在仓库根目录的 MATLAB 命令窗口运行：

```matlab
addpath('matlab/resonance_qfit');
run_qfit
```

在 `run_qfit.m` 顶部选择幅值/功率/dB/复数和线形。支持 Lorentzian、Fano、Gaussian、Voigt、pseudo-Voigt。`fit_resonance_q.m` 可直接处理数组，`batch_qfit.m` 可批处理。

不依赖额外拟合工具箱；使用基础 MATLAB `fminsearch`，多初值优化。演示：`demo_qfit`；数值和读写测试：`test_qfit`。

Fano 的拟合极点 Q 与可见峰/谷宽度 Q 可能不同。完整用法、公式、输出定义和限制见 [说明文档](../../docs/resonance_qfit.md)。
