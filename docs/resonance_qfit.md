# 共振线形拟合与 Q 值计算

本工具面向 CST 导出的频率响应以及实验频谱，提供 Python 和 MATLAB 两套实现。支持 **Lorentzian、Fano、Gaussian、Voigt、pseudo-Voigt**，支持峰和谷、常数或线性局部背景、模型比较和批处理。不做平滑，不替换尖峰，不从截图自动生成“真值”。

## 先选对数据类型

| 导出的量 | `InputKind` / `--kind` | 程序使用的功率谱 |
|---|---|---|
| CST Magnitude，线性幅值 `abs(S21)` | `amplitude` | `abs(S21)^2` |
| 功率透射率或反射率 | `power` | 原值 |
| S 参数 dB：`20*log10(abs(S21))` | `db` | `10^(dB/10)` |
| 复数 S21，实部和虚部 | `complex` | `abs(S21)^2`；相位未参与拟合 |

功率减半等于幅值乘 `1/sqrt(2)`，不是幅值减半。例如幅值峰值为 0.94、背景可忽略时，3 dB 水平约为 0.665。

频率列必须是真正的频率，单位保持一致即可；Hz、GHz、THz 得到相同的无量纲 Q。不要直接输入波长或时间列。Python 默认列号 0、1；MATLAB 默认列号 1、2。

一个文件必须对应一个几何/一次扫描的一条谱。TXT/CSV/TSV 支持 UTF-8 文本表头、注释和 CST 分隔行，十进制小数点为 `.`；不支持十进制逗号、Excel 和带引号的复杂 CSV。重复频率会报错，防止将不同 CST 参数块平均成一条假谱。

## 为什么返回三种 Q

在选定的局部窗口内拟合 **功率谱**，背景 `B(f)` 可为常数或线性函数。

Lorentzian：

```text
P(f) = B(f) + A / (1 + epsilon^2)
epsilon = 2*(f-f0)/Gamma
```

Fano 采用数值更稳定的等价表达式：

```text
P(f) = B(f) + (C + D*epsilon)/(1 + epsilon^2)
epsilon = 2*(f-f0)/Gamma
```

这对应常见的 `B' + S*(q+epsilon)^2/(1+epsilon^2)`，其中常数项并入背景。输出的 `fano_q_canonical_positive_scale` 按 `S>=0` 的约定选支；在 Lorentz 极限为 `null`（MATLAB 为 `[]`）。有损或复杂背景下，它只是所选强度模型的有效 q，不能单凭该值确定物理机制。

| 输出字段 | 定义 | 适用范围 |
|---|---|---|
| `q_resonance` | `f0/Gamma` | Lorentzian、单共振 Fano 的拟合极点线宽 Q；Fano 的 f0 一般不等于可见峰位 |
| `q_spectral` | `feature_frequency/feature_fwhm` | 从拟合曲线减去背景，取选定峰/谷相对背景的半高/半深宽 |
| `q_3db` | 功率峰频率/绝对半功率带宽 | 仅正峰；用完整拟合功率曲线的峰值一半，要求窗口内两边都有交点 |

`q_3db` 的交点水平是 **完整峰值的一半**，`q_spectral` 的交点水平是 **相对局部背景的对比度的一半**。非零背景时两者会不同。缺少两侧交点就返回空值，并提示扩大窗口或检查定义；程序不会把窗口端点冒充交点。

Fano 的可见峰/谷半宽通常不等于 Gamma。示例中 f0=0.33 THz、Gamma=0.006 THz，所以 `q_resonance=55`；可见峰位和 `q_spectral` 可以明显不同。这不代表程序出错，比较数据集时必须固定 Q 定义。

Gaussian 的宽度是功率高斯 FWHM；Voigt 是 Gaussian 与 Lorentzian 的卷积；pseudo-Voigt 是共用 FWHM 的两种线形加权和。后三者返回 `q_spectral`，`q_resonance` 留空。Voigt 额外返回 Gaussian/Lorentzian 分量宽度和 `q_lorentz_component`，只有在仪器高斯展宽等假设成立、分量可辨识时才有相应物理含义。光谱拟合不能自动分离辐射 Q、吸收 Q 或内禀 Q。

## Python 使用

Python 3.10+，依赖 NumPy、SciPy、Matplotlib。在仓库根目录运行：

```bash
python -m pip install -r python/resonance_qfit/requirements.txt
python python/resonance_qfit/qfit.py "D:/spectra/S21.txt" --kind amplitude --model fano --feature peak --window 0.30 0.36 --output results/sample1
```

**将窗口改成实际目标共振附近的范围。** 相邻存在多个模式时，先分窗；不要把整段多模态频谱强行拟合为一个 Fano。保留峰两侧足够背景，同时避免把下一个模式放进窗口。

输出：`sample1.json` 参数/Q/诊断、`sample1.csv` 原始功率/拟合/背景/残差、`sample1.png` 拟合与残差图。`--no-plot` 可跳过画图。不同输入数据用不同输出前缀；命令行会阻止直接覆盖输入文件。

其他用法：

```bash
# 负峰/透射谷
python python/resonance_qfit/qfit.py dip.csv --kind power --model lorentzian --feature dip --window 0.30 0.36 --output results/dip

# 将五种线形用同一窗口的 AICc 排序
python python/resonance_qfit/qfit.py S21.txt --kind amplitude --model auto --window 0.30 0.36 --output results/compare

# 频率、S21实部、S21虚部三列；Python列号从0开始
python python/resonance_qfit/qfit.py complex.txt --kind complex --imag-col 2 --model fano --window 0.30 0.36 --output results/complex

# 一个文件一条谱，批量提取同一模式；不同模式应分别设置窗口
python python/resonance_qfit/batch_qfit.py "D:/spectra" --pattern "*.txt" --kind amplitude --model fano --feature peak --window 0.30 0.36 --output-dir results/batch --no-plots

# 合成数据演示与测试
python python/resonance_qfit/demo.py
python -m unittest discover -s python/resonance_qfit/tests -v
```

批处理生成 `summary.csv`；单文件错误会记入 `error` 列，其他文件继续运行。`fit_completed` 只表示拟合已完成，不能代替诊断和科学审核。

作为函数调用：

```python
import sys
sys.path.insert(0, "python/resonance_qfit")
from qfit import load_spectrum, fit_resonance

f, amplitude = load_spectrum("S21.txt")
r = fit_resonance(f, amplitude, kind="amplitude", model="fano",
                  feature="peak", window=(0.30, 0.36))
print(r.f0, r.q_resonance, r.q_spectral, r.warnings)
```

## MATLAB 使用

使用基础 MATLAB 的 `fminsearch`，不依赖 Optimization Toolbox 或 Curve Fitting Toolbox。Voigt 用复误差函数的有理展开计算，方法见 Weideman 1994。界面/导出功能建议 MATLAB R2020a+。

```matlab
addpath('matlab/resonance_qfit');
run_qfit                  % 选择文本谱，输入拟合频率上下限
```

运行前在 `run_qfit.m` 顶部设置 `cfg.inputKind`、`cfg.model` 和 `cfg.feature`。默认 `amplitude / fano / peak`。输出在输入文件夹下的 `qfit_results/`，含时间戳。

函数方式：

```matlab
[f, amp] = read_cst_spectrum('D:/spectra/S21.txt', 1, 2);
r = fit_resonance_q(f, amp, 'Model','fano', 'InputKind','amplitude', ...
    'Feature','peak', 'Window',[0.30 0.36], 'FrequencyUnit','THz');
disp(r.q_resonance);
disp(r.q_spectral);
disp(r.warnings);
save_qfit_result(r, 'results/sample1_matlab');
% 无图形环境可仅导出 JSON/CSV：save_qfit_result(r,'results/sample1_matlab',false);

% 复数：第2列实部，第3列虚部
[f, s21] = read_cst_spectrum('complex.txt',1,2,3);
r = fit_resonance_q(f,s21,'InputKind','complex','Model','fano','Window',[.30 .36]);

% 批处理（默认频率/信号为第1/2列）
summary = batch_qfit('D:/spectra','*.txt','D:/fits', ...
    'Model','fano','InputKind','amplitude','Feature','peak','Window',[.30 .36]);

demo_qfit;
test_qfit;
```

MATLAB 用 `[]` 表示不可用数值，JSON 中为 `[]`；Python 为 `None`，JSON 中为 `null`。

## 使用前必须检查的结果

- 实际求解频点能否解析线宽。少于约 8 个采样间隔/相关线宽会提示欠采样；插值导出更密的点无法增加真实分辨率。
- 左右半高/半深交点是否都位于窗口内。截断峰不应给出“高 Q 已证实”的判断。
- 残差是否有连续结构或额外峰谷。均方误差小或 R² 很高也可能得到错误的线宽；窗口过窄时背景与宽度还可能相互补偿。
- 换一个合理窗口、增加初值、加密仿真网格后，f0 和 Q 是否稳定。
- Fano/Gaussian/Voigt 的 AICc 排名只是候选模型比较，不证明 Fano 机制、BIC、准 BIC 或仪器展宽的存在。

本版本针对**单个可分离共振**。复杂 EIT、多个重叠极点、强起伏玻璃背景，需要相应的耦合模/多共振模型，不能因为软件给了一个数字就认定是该模式的物理 Q。当前没有置信区间、全局最优保证或自动异常谱剔除。保留原谱和求解设置，避免将未收敛数据当训练标签。

## 方法来源

本次验证（2026-10-26）：Python 3.13、NumPy 2.1.3、SciPy 1.15.3、Matplotlib 3.10.0 的 14 项测试通过；MATLAB R2025a 的五类线形、峰/谷、单位与 dB/复数转换、自动比较、欠采样/截窗提示、CST 风格文本、JSON/CSV/PNG 和批处理检查通过。MATLAB Voigt 测试以独立数值卷积积分作参照。测试输入均为合成数据；尚未用你的实际 CST 导出文件做误差评估，交互选文件对话框未做自动化验证。

- Fano 功率线形及 `epsilon=2*(f-f0)/Gamma`、`Q=f0/Gamma` 定义：[Controlling Fano resonances in multilayer dielectric gratings towards optical bistable devices](https://doi.org/10.1038/s41598-018-34787-9)。
- Lorentzian/Gaussian/Voigt/pseudo-Voigt 定义：[lmfit 官方模型文档](https://lmfit.github.io/lmfit-py/builtin_models.html)。本代码没有依赖 lmfit。
- Python Voigt：[SciPy `voigt_profile`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.voigt_profile.html)。
- MATLAB 复误差函数算法：J. A. C. Weideman, *Computation of the Complex Error Function*, SIAM J. Numer. Anal. 31, 1497–1518 (1994), [DOI](https://doi.org/10.1137/0731077)。系数在代码中生成。
