# Python 共振 Q 值工具

第一次使用请看 [逐步指南](../../docs/START_HERE.md)。仓库根目录的 `fit_q.py` 是较短的命令入口。

支持 Lorentzian、Fano、Gaussian、Voigt、pseudo-Voigt；幅值/功率/dB/复数输入，峰/谷拟合、自动模型比较和批处理。

在仓库根目录运行：

```bash
python -m pip install -r python/resonance_qfit/requirements.txt
python python/resonance_qfit/qfit.py S21.txt --kind amplitude --model fano --feature peak --window 0.30 0.36 --output results/sample1
```

窗口请改为实际目标共振附近。Fano 的 `q_resonance` 与可见峰/谷的 `q_spectral` 一般不同。高斯/Voigt 线形的谱宽 Q 不能直接当作腔衰减 Q。

完整用法、公式、输出定义和限制见 [说明文档](../../docs/resonance_qfit.md)。
