# MATLAB 共振 Q 值工具

交互计算优先使用根目录的 **[qcalc_gui.m](../../qcalc_gui.m)**，只下载该文件即可运行，不需要本目录。窗口中可选文件、频率区间、线形、背景和峰/谷，并导出结果。操作见 [使用指南](../../docs/START_HERE.md)。

本目录保留不同用途的函数接口、批处理和验证：

- `fit_resonance_q.m`：直接拟合数组，支持 Fano、Lorentzian、Gaussian、Voigt、pseudo-Voigt。
- `qfit_file.m`：以文件为输入调用拟合并保存结果。
- `batch_qfit.m`：批量处理单组频谱文件。
- `read_cst_spectrum.m`、`save_qfit_result.m`：读写辅助。
- `demo_qfit.m`、`test_qfit.m`：合成数据演示及数值检查。

```matlab
addpath('matlab/resonance_qfit');
[f, amp] = read_cst_spectrum('S21.txt');
r = fit_resonance_q(f, amp, 'Model','fano', 'Window',[.46 .49]);
disp(r.q_resonance);
```

示例窗口应改成你的实际共振区间。使用基础 MATLAB `fminsearch`，无需额外拟合工具箱。Fano 极点 Q 与可见峰宽 Q 不同；完整定义见 [说明文档](../../docs/resonance_qfit.md)。
