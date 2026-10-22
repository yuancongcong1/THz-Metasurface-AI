clear; clc; close all;

%% ================= 用户设置 =================

timeCol = 1;                 % 时间列
refCol  = 3;                 % 参考信号列
samCol  = 7;                 % 样品信号列

inputTimeUnit = 'auto';      % 's' 或 'ps'；auto 识别表头/THz采样间隔
inputFieldUnit = 'auto';     % 'V/m'、'V'、'a.u.' 等；auto 读取模拟数据表头
inputFile = '';              % 留空弹窗选文件；也可填写文件完整路径

fftStart = 0e-12;               % FFT起始时间，s；[]表示不限
fftEnd   = 120e-12;          % FFT结束时间，s；[]表示不限

zeroNum = 0;              % 末尾补零点数，0表示不补零
removeMean = false;          % 是否去除截取信号的平均值

useGaussianWindow = true;   % 高斯窗开关：true 开启，false 关闭
gaussianCenter = [54.7e-12,56.8e-12];        % s；[] 分别自动寻中心；[参考, 样品] 分别指定；标量指定共同中心
gaussianSigma = [];         % s；[] 自动设宽；[参考, 样品] 分别指定；标量指定共同标准差
gaussianWidthMode = 'matched'; % 自动设宽：matched 同宽异中心；independent 独立宽度
gaussianWidthFactor = 3.5;  % 自动标准差/各自主脉冲半宽；越大保留的时域范围越宽
% 自动寻窗只在 fftStart 至 fftEnd 的截取区间内进行。
% 中心按时域主脉冲位置确定，不按频谱峰频率确定。
% 比较透射谱时默认同宽异中心；独立宽度会引入不同的频谱展宽。
% gaussianWidthMode 只控制自动宽度；填写 gaussianSigma 时以手动值为准。

freqRange = [0.5e12, 1.2e12]; % 频谱显示范围，Hz
refThreshold = 0.01;         % 参考频谱低于最大幅度1%的位置不画透射率

%% ================= 读取文件 =================

if isempty(inputFile)
    [file, folder] = uigetfile( ...
        {'*.xlsx;*.xls;*.csv;*.txt;*.dat', '时域数据文件'}, ...
        '选择数据文件：时间、参考、样品');
    if isequal(file, 0)
        return;
    end
    filename = fullfile(folder, file);
else
    filename = inputFile;
end

assert(isfile(filename), '找不到数据文件：%s', filename);
[~, ~, ext] = fileparts(filename);

if ismember(lower(ext), {'.xlsx', '.xls'})
    data = readmatrix(filename, 'Sheet', 1);
    header = readcell(filename, 'Sheet', 1, 'Range', 'A1:C1');
else
    data = readmatrix(filename, 'FileType', 'text');
    header = {};
end

assert(size(data, 2) >= max([timeCol, refCol, samCol]), ...
    '文件列数不足，请检查列号。');

data = data(:, [timeCol, refCol, samCol]);

% 去掉表头或首尾空行产生的无效行
valid = all(isfinite(data), 2);
first = find(valid, 1, 'first');
last  = find(valid, 1, 'last');

assert(~isempty(first), '没有读取到有效数据。');
data = data(first:last, :);

assert(all(isfinite(data(:))), ...
    '数据中间存在空值或非数字，请检查文件。');

tInput = data(:, 1);
E_ref = data(:, 2);
E_sam = data(:, 3);

%% ================= 检查采样 =================

assert(numel(tInput) >= 3, '数据点不足。');
assert(all(diff(tInput) > 0), '时间必须严格递增。');

% 文件边界处确定时间单位；之后所有时间量均为 s
resolvedTimeUnit = lower(inputTimeUnit);
if strcmp(resolvedTimeUnit, 'auto')
    if numel(header) >= timeCol && ...
            strcmpi(string(header{timeCol}), 'Time_s')
        resolvedTimeUnit = 's';
    else
        dtInput = median(diff(tInput));
        if dtInput > 0 && dtInput <= 1e-9
            resolvedTimeUnit = 's';
        elseif dtInput >= 1e-3 && dtInput <= 10
            resolvedTimeUnit = 'ps';
            warning('无时间单位表头：按 THz 数据的采样间隔推定时间列为 ps；如不符，请设置 inputTimeUnit。');
        else
            error('无法可靠判断时间单位；请将 inputTimeUnit 设为 ''s'' 或 ''ps''。');
        end
    end
end

switch resolvedTimeUnit
    case 's'
        t = tInput;
    case 'ps'
        t = tInput * 1e-12;
    otherwise
        error('inputTimeUnit 必须为 ''auto''、''s'' 或 ''ps''。');
end

if strcmpi(inputFieldUnit, 'auto')
    if numel(header) >= max([refCol, samCol]) && ...
            strcmpi(string(header{refCol}), 'E_ref_V_per_m') && ...
            strcmpi(string(header{samCol}), 'E_sam_V_per_m')
        fieldUnit = 'V/m';
    else
        fieldUnit = 'a.u.'; % 无标定信息时不假设绝对电场单位
    end
else
    fieldUnit = inputFieldUnit;
end

fprintf('输入时间单位：%s；内部时间单位：s；信号单位：%s\n', ...
    resolvedTimeUnit, fieldUnit);

dt = mean(diff(t));

assert(max(abs(diff(t)-dt)) <= 1e-3*dt, ...
    '时间采样不均匀，请检查时间列。');

assert(isscalar(zeroNum) && zeroNum >= 0 && ...
    zeroNum == fix(zeroNum), '补零数必须是非负整数。');

%% ================= 截取信号 =================

idx = true(size(t));
timeTolerance = dt*1e-6;     % 避免浮点舍入漏掉窗口端点

if ~isempty(fftStart)
    idx = idx & t >= fftStart-timeTolerance;
end

if ~isempty(fftEnd)
    idx = idx & t <= fftEnd+timeTolerance;
end

assert(nnz(idx) >= 3, '截取范围内数据点不足。');

t_cut = t(idx);
ref_cut = E_ref(idx);
sam_cut = E_sam(idx);

if removeMean
    ref_cut = ref_cut - mean(ref_cut);
    sam_cut = sam_cut - mean(sam_cut);
end

%% ================= 可选高斯窗 =================

assert(isscalar(useGaussianWindow) && ...
    (islogical(useGaussianWindow) || isnumeric(useGaussianWindow)) && ...
    ismember(useGaussianWindow, [0, 1]), ...
    'useGaussianWindow 必须为 true 或 false。');

fftWindows = ones(numel(t_cut), 2); % 两列分别为参考窗、样品窗，无量纲
windowCenters = nan(1, 2);         % 两列顺序均为 [参考, 样品]，s
windowSigmas = nan(1, 2);          % s
windowDescription = 'No window';

if useGaussianWindow
    configuredCenters = expandGaussianSetting(gaussianCenter, 'gaussianCenter');
    configuredSigmas = expandGaussianSetting(gaussianSigma, 'gaussianSigma');
    assert(all(isnan(configuredCenters) | ...
        (configuredCenters >= t_cut(1) & configuredCenters <= t_cut(end))), ...
        'gaussianCenter 必须位于截取区间内，单位 s。');
    assert(all(isnan(configuredSigmas) | configuredSigmas > 0), ...
        'gaussianSigma 必须为正数，单位 s。');
    assert(any(strcmpi(gaussianWidthMode, {'matched', 'independent'})), ...
        'gaussianWidthMode 必须为 ''matched'' 或 ''independent''。');
    assert(isnumeric(gaussianWidthFactor) && isscalar(gaussianWidthFactor) && ...
        isreal(gaussianWidthFactor) && isfinite(gaussianWidthFactor) && ...
        gaussianWidthFactor >= 1, 'gaussianWidthFactor 必须为不小于 1 的有限实数。');

    autoPulseRanges = nan(2, 2);
    pulseDetected = false(2, 1);
    cutSignals = {ref_cut, sam_cut};
    windowNames = {'Reference', 'Sample'};

    for k = 1:2
        if isnan(configuredCenters(k)) || isnan(configuredSigmas(k))
            [autoPulseRanges(k,:), pulseDetected(k)] = detectMainPulse(t_cut, cutSignals{k});
            if ~pulseDetected(k)
                warning('%s 未检测到明显主脉冲，自动参数回退到截取区间中心和时长/6。', ...
                    windowNames{k});
            elseif autoPulseRanges(k,1) <= t_cut(1)+timeTolerance || ...
                    autoPulseRanges(k,2) >= t_cut(end)-timeTolerance
                warning('%s 主脉冲接近截取边界，请检查 FFT 截取区间。', windowNames{k});
            end
        end

        if isnan(configuredCenters(k))
            if pulseDetected(k)
                windowCenters(k) = mean(autoPulseRanges(k,:));
            else
                windowCenters(k) = t_cut(1) + (t_cut(end)-t_cut(1))/2;
            end
        else
            windowCenters(k) = configuredCenters(k);
        end

        if isnan(configuredSigmas(k))
            if pulseDetected(k)
                pulseHalfWidth = max(abs(autoPulseRanges(k,:)-windowCenters(k)));
                windowSigmas(k) = gaussianWidthFactor * max(pulseHalfWidth, dt);
            else
                windowSigmas(k) = (t_cut(end)-t_cut(1))/6;
            end
        else
            windowSigmas(k) = configuredSigmas(k);
        end
    end

    % 共同标准差取两条自动估计值中较大者；中心保持各自独立。
    if isempty(gaussianSigma) && strcmpi(gaussianWidthMode, 'matched')
        windowSigmas(:) = max(windowSigmas);
    end

    % 直接按时间计算，无需 Signal Processing Toolbox
    for k = 1:2
        fftWindows(:,k) = exp(-0.5*((t_cut-windowCenters(k))/windowSigmas(k)).^2);
        assert(max(fftWindows(:,k)) > 1e-12, ...
            '%s 高斯窗过窄，采样点上的窗值接近零。', windowNames{k});
        fprintf('%s 高斯窗：中心 %.6g s，标准差 %.6g s\n', ...
            windowNames{k}, windowCenters(k), windowSigmas(k));
    end
    windowDescription = 'Separate Gaussian windows';
end

% 在原时间轴上分别加窗，再补零；不移动脉冲、不重置各自的时间原点。
fftWindowRef = fftWindows(:,1);
fftWindowSam = fftWindows(:,2);
ref_fft = ref_cut .* fftWindowRef;
sam_fft = sam_cut .* fftWindowSam;
% 不对脉冲数据自动做相干增益补偿：加窗后的幅值代表加权信号。
% 即使同宽异中心，加窗仍可能影响不同形状脉冲的幅值比，应比较开关前后结果。

%% ================= 时域图 =================

timeFigure = figure('Color', 'w');

signals = {E_ref, E_sam};
names = {'Reference', 'Sample'};
fftSignals = {ref_fft, sam_fft};

for k = 1:2
    subplot(2, 1, k);
    hold on;

    E = signals{k};

    margin = 0.12 * (max(E)-min(E));
    if margin == 0
        margin = max(1, abs(E(1))) * 0.1;
    end

    yLow = min(E)-margin;
    yHigh = max(E)+margin;

    % 绿色透明底色：完整原始信号
    hOriginal = patch( ...
        [t(1), t(end), t(end), t(1)], ...
        [yLow, yLow, yHigh, yHigh], ...
        [0.2, 0.8, 0.3], ...
        'FaceAlpha', 0.12, ...
        'EdgeColor', [0, 0.55, 0], ...
        'LineWidth', 1.2);

    % 红色透明底色：参与FFT的部分
    hFFT = patch( ...
        [t_cut(1), t_cut(end), t_cut(end), t_cut(1)], ...
        [yLow, yLow, yHigh, yHigh], ...
        [1, 0.2, 0.2], ...
        'FaceAlpha', 0.22, ...
        'EdgeColor', [0.85, 0, 0], ...
        'LineWidth', 1.2);

    hSignal = plot(t, E, 'b', 'LineWidth', 1);

    xlabel('Time (s)');
    ylabel(['Signal amplitude (' fieldUnit ')']);
    title(names{k});

    if useGaussianWindow
        hWindowed = plot(t_cut, fftSignals{k}, 'm--', 'LineWidth', 1.2);
        % 高斯窗是无量纲权重，使用独立右轴，避免与电场单位混淆。
        yyaxis right;
        hGaussian = plot(t_cut, fftWindows(:,k), 'k-.', 'LineWidth', 1.6);
        ylabel('Gaussian window weight (dimensionless)');
        ylim([0, 1.05]);
        ax = gca;
        ax.YAxis(2).Color = [0.15, 0.15, 0.15];
        yyaxis left;
        ax.YAxis(1).Color = [0, 0.25, 0.7];
        legend([hSignal, hOriginal, hFFT, hWindowed, hGaussian], ...
            {'Signal', 'Original range', 'FFT range', ...
             'Windowed FFT input', 'Gaussian window w(t)'}, ...
            'Location', 'best');
        title(sprintf('%s: center = %.4g s, sigma = %.4g s', ...
            names{k}, windowCenters(k), windowSigmas(k)));
    else
        legend([hSignal, hOriginal, hFFT], ...
            {'Signal', 'Original range', 'FFT range'}, ...
            'Location', 'best');
    end

    xlim([t(1), t(end)]);
    ylim([yLow, yHigh]);
    grid on;
end

if useGaussianWindow
    sgtitle('Separate Gaussian windows for reference and sample');
end

%% ================= 补零和FFT =================

N = numel(t_cut);
Nfft = N + zeroNum;

% 自动在信号末尾补零
F_ref = fft(ref_fft, Nfft);
F_sam = fft(sam_fft, Nfft);

nPositive = floor(Nfft/2)+1;

% dt 的单位为 s，故频率轴的单位为 Hz
f = (0:nPositive-1)' / (Nfft*dt);

% 单边幅度谱：按补零前点数归一化；开启高斯窗后为加权信号的幅度谱
amp_ref = abs(F_ref(1:nPositive)) / N;
amp_sam = abs(F_sam(1:nPositive)) / N;

if mod(Nfft, 2) == 0
    amp_ref(2:end-1) = 2*amp_ref(2:end-1);
    amp_sam(2:end-1) = 2*amp_sam(2:end-1);
else
    amp_ref(2:end) = 2*amp_ref(2:end);
    amp_sam(2:end) = 2*amp_sam(2:end);
end

%% ================= FFT幅度谱 =================

figure('Color', 'w');

plot(f, amp_ref, 'b', 'LineWidth', 1.2);
hold on;
plot(f, amp_sam, 'r', 'LineWidth', 1.2);

xlabel('Frequency (Hz)');
ylabel(['FFT amplitude (' fieldUnit ')']);
legend('Reference', 'Sample');
title(['FFT spectra - ' windowDescription]);
xlim(freqRange);
grid on;

%% ================= 幅度透射谱 =================

T = nan(size(f));

valid = amp_ref > refThreshold*max(amp_ref);
T(valid) = amp_sam(valid) ./ amp_ref(valid);

figure('Color', 'w');

plot(f, T, 'k', 'LineWidth', 1.2);

xlabel('Frequency (Hz)');
ylabel('Amplitude transmission (dimensionless)');
title(['Sample / Reference - ' windowDescription]);
xlim(freqRange);
grid on;

%% ================= 高斯窗参数展开 =================
function values = expandGaussianSetting(setting, settingName)
% []：两条信号均自动；标量：共同手动值；二元向量：[参考, 样品]。
if isempty(setting)
    values = [NaN, NaN];
    return;
end
assert(isnumeric(setting) && isreal(setting) && isvector(setting) && ...
    ismember(numel(setting), [1, 2]) && all(isfinite(setting(:))), ...
    '%s 必须为 []、有限实数标量或 [参考, 样品] 二元向量。', settingName);
values = reshape(setting, 1, []);
if isscalar(values)
    values = [values, values];
end
end

%% ================= 自动主脉冲检测（无需额外工具箱） =================
function [pulseRange, found] = detectMainPulse(tSignal, signal)
% 用解析信号包络定位最强主脉冲，不将后续独立回波并入主脉冲范围。
% 检测用的去基线、归一化和包络平滑不改变传入 FFT 的信号。
pulseRange = [NaN, NaN];
found = false;
x = signal(:) - median(signal);
scale = max(abs(x));
if scale == 0
    return;
end
x = x/scale;
sampleCount = numel(x);

% FFT 构造解析信号；补零降低首尾周期连接对包络的影响。
analyticLength = 2^nextpow2(2*sampleCount);
analyticWeights = zeros(analyticLength, 1);
analyticWeights(1) = 1;
analyticWeights(2:analyticLength/2) = 2;
analyticWeights(analyticLength/2+1) = 1;
analyticSignal = ifft(fft(x, analyticLength).*analyticWeights);
envelope = abs(analyticSignal(1:sampleCount));
smoothPoints = max(3, 2*floor((sampleCount-1)/400)+1);
envelope = movmean(envelope, smoothPoints);

% 用较低的一半包络值估计背景，降低独立回波对背景估计的影响。
sortedEnvelope = sort(envelope);
background = sortedEnvelope(1:ceil(sampleCount/2));
noiseFloor = median(background);
noiseSpread = 1.4826*median(abs(background-noiseFloor));
[peakValue, peakIndex] = max(envelope);
if peakValue-noiseFloor <= max(6*noiseSpread, 0.1*peakValue)
    return; % 连续振荡或没有明显脉冲时，不强行指定一个局域主脉冲。
end
threshold = min(0.8*peakValue, max(0.15*peakValue, noiseFloor+3*noiseSpread));

left = peakIndex;
while left > 1 && envelope(left-1) >= threshold
    left = left-1;
end
right = peakIndex;
while right < sampleCount && envelope(right+1) >= threshold
    right = right+1;
end
% 轻微外扩，补偿包络平滑带来的边界不确定性。
padding = floor(smoothPoints/2);
left = max(1, left-padding);
right = min(sampleCount, right+padding);
if right-left < 2 || right-left >= 0.9*(sampleCount-1)
    return;
end
pulseRange = [tSignal(left), tSignal(right)];
found = true;
end


