# 第一次使用这个仓库

先做最容易的一件事：用 `p260r110.txt` 算 Q。完整文件夹需要保留，不要只下载一个 `.m` 文件。

## 一 先拿到代码

仓库现为公开仓库。浏览器打开：

https://github.com/yuancongcong1/THz-Metasurface-AI

点击绿色 `Code` → `Download ZIP`，再右键解压 ZIP。建议把解压后包含 `START_HERE.m` 的文件夹放到：

```text
E:\YCC\THz-Metasurface-AI\
```

位置可以不同；入口会自动寻找自身文件夹。不要在压缩包预览里直接运行文件。

“离线代码包”就是同样的源码 ZIP。把它复制到另一台服务器后解压，运行时不需要连接 GitHub；服务器仍需安装 MATLAB，或者 Python 及相应依赖。ZIP 中没有个人实测/CST 数据、MATLAB、CST 或许可证。

## 二 用你的 p260r110 数据运行 MATLAB

准备好数据文件：

```text
E:\YCC\20261026\p260r110.txt
```

1. 打开 MATLAB。
2. 点击“打开”，选择代码文件夹根目录中的 `EXAMPLE_P260R110.m`。
3. 点击编辑器顶部的绿色“运行”按钮。如果 MATLAB 问是否更改文件夹或添加路径，选择更改文件夹即可。
4. 程序会读取上面的 E 盘路径；如果该电脑没有这个路径，会弹窗让你选择 `p260r110.txt`。
5. 等命令窗口打印结果。程序会自动设置 **Magnitude 幅值、Fano、透射峰、THz、拟合窗口 0.466–0.471 THz**，并弹出拟合/残差图。

这个示例窗口只适用于这份示例文件。之前文档中的 0.30–0.36 THz 是格式演示，对频段 0.45–0.50 THz 的文件会导致窗口内没有数据。

输出位于数据旁边新建的文件夹：

```text
E:\YCC\20261026\qfit_results\
  p260r110_日期时间.json  参数、Q值、提示
  p260r110_日期时间.csv   输入功率、拟合功率、背景、残差
  p260r110_日期时间.png   拟合图和残差图
```

`result.q_resonance` 是 Fano 的拟合极点线宽 Q；`result.q_spectral` 是减去拟合背景后可见峰的半高宽 Q；`result.q_3db` 是完整功率峰下降一半的带宽 Q。三者定义不同，不应混用。参考文献用哪种定义，就选相同的列。

拟合得到的 Q 只描述导出曲线。导出点很多不代表 CST 在每一点都实际求解，也不能证明很尖的峰已完成网格/频率收敛验证。Fano 拟合不能识别或排除 BIC/准 BIC。

## 三 换成另一份谱时怎样操作

打开根目录 `START_QFIT.m`，点击绿色“运行”，依次操作：

1. 选择 TXT/CSV/TSV/DAT 文件。
2. 选择数据类型。CST 图标题为 `S-Parameters [Magnitude]` 时选第一项 **Magnitude**；已经平方的透射率选 Power；dB 曲线选第三项。
3. 确认列号。通常第 1 列频率、第 2 列信号，单位 THz。
4. 频谱图出现后，可以用缩放功能放大一个目标共振。点图下方“放大完成，开始点击左右边界”按钮。
5. 在共振左侧、右侧各点击一次。对 Fano，窗口应同时包含相邻峰和谷，并留出两侧背景。
6. 在弹窗中核对上下限，再选 Fano、Lorentzian 或 Auto。
7. 选择“峰”或“谷”，等待输出。

多个远离的共振分开拟合。不要用单个 Fano 拟合所有玻璃起伏和多个模式。

## 四 Python 的完整命令

在 Windows 代码文件夹空白处右键“在终端中打开”，依次执行：

```powershell
python -m pip install -r python/resonance_qfit/requirements.txt
python fit_q.py "E:\YCC\20261026\p260r110.txt" --kind amplitude --model fano --feature peak --window 0.466 0.471 --output "E:\YCC\20261026\qfit_results\p260r110_python"
```

第一行安装一次即可；第二行运行你的数据。`--output` 后面是输出前缀，不需要写 `.png` 后缀。相同前缀重跑会更新对应结果，请为要保留的版本换一个前缀。

若电脑提示找不到 `python`，先安装 Python/Anaconda，或打开已安装环境的终端。若提示 `No module named scipy`，在实际运行的那个 Python 环境执行第一行安装命令。

## 五 在其他服务器上用

现在可以匿名下载公开代码。服务器有网络和 Git 时：

```bash
git clone https://github.com/yuancongcong1/THz-Metasurface-AI.git
cd THz-Metasurface-AI
python3 -m pip install -r python/resonance_qfit/requirements.txt
python3 fit_q.py "/home/你的账号/data/p260r110.txt" --kind amplitude --model fano --feature peak --window 0.466 0.471 --output results/p260r110
```

最后一行的输入路径要改为数据在服务器上的真实位置。Windows 的 `E:\...` 路径不会自动出现在 Linux 服务器上；需要把 TXT 数据也复制过去。公开仓库只解决代码访问权限，不解决服务器断网、DNS/代理或 Python 依赖缺失。

没有网络时，将代码 ZIP 和你自己的 TXT 文件分别上传到服务器再解压。MATLAB 路线需要服务器有 MATLAB；Python 路线需要服务器已有 NumPy、SciPy、Matplotlib，或者另外准备匹配其系统/Python版本的离线依赖包。

## 六 仓库里其他脚本分别做什么

| 想做的事 | 文件 | 操作前需要设置 | 得到什么 |
|---|---|---|---|
| 打开任务菜单 | 根目录 `START_HERE.m` | 无 | 选择任务的入口 |
| 拟合一条 CST 谱、算 Q | `START_QFIT.m` / `EXAMPLE_P260R110.m` | 数据类型、一个共振的窗口 | Q、拟合图、残差和参数 |
| 一批 TXT 谱算 Q | `matlab/resonance_qfit/batch_qfit.m` 或 Python 的 `batch_qfit.py` | 输入/输出目录、同一模式窗口 | 每条谱的结果和 summary.csv |
| 生成模拟 THz 时域信号 | `matlab/tds/data_sin_shengcheng.m` | 可先保留默认示例 | 人工生成的参考/样品信号与 THz_simulated.xlsx |
| 处理实际 THz-TDS 时域数据 | `matlab/tds/shujuchuli.m` | 时间/参考/样品列号、单位、FFT区间、是否加窗 | 时域、FFT幅度谱和幅度透射谱 |
| 生成 CST 参数表 | `matlab/cst/generate_parameters.m` | 你的结构参数名、范围、几何/加工约束 | parameters.tsv；本步骤不运行 CST |
| 批量修改参数并运行 CST | `matlab/cst/run_cst_sweep.m` | CST工程、参数表、求解器和结果树路径 | 每组几何的 S21/S11 文件及日志 |
| 画结构示意图 | `visualization/render_metasurface.py` | Blender、示意图参数和输出目录 | PNG、可编辑 .blend、.glb；不是电磁仿真结果 |

TDS 示例生成脚本会在当前工作文件夹输出 Excel；从新菜单生成时输出到仓库 `results/tds_demo/`。处理该模拟三列文件时，`shujuchuli.m` 顶部设：

```matlab
timeCol=1; refCol=2; samCol=3;
inputTimeUnit='s';
fftStart=0; fftEnd=20e-12;
useGaussianWindow=false;
```

处理实测文件时按实测列号和时域范围设置。不能把默认第 1/3/7 列用于所有文件。

CST 两个脚本当前仍是 **AuSRR 模板**。它们默认的 `p,r,w,d,l1,l2` 参数及 20–30 μm 线宽范围，不适用于你的三开口环槽和大于 30 μm 的加工约束，必须先改。CST 批处理默认结果树使用 mode 2，而你的截图使用 mode 1；也要按实际结果树修改。若导出的是 Magnitude，应使用 `linear`，不能再当 dB 转换。Windows 电脑上需要同时安装 MATLAB 和 CST。

Blender 示意图从仓库根目录执行：

```bash
blender -b --factory-startup --python visualization/render_metasurface.py -- --output-dir outputs --style paper --view plan
```

目前仓库尚未提供神经网络正向/逆向训练脚本；参数生成、CST 数据采集、Q 提取和绘图是为后续训练准备的工具。

更多 Q 定义与拟合限制见 [详细文档](resonance_qfit.md)。
