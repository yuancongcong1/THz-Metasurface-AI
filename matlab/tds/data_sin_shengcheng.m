clear; clc; close all;

%% 1. 参数设置
ps = 1e-12;                   % 1 ps = 1e-12 s
THz = 1e12;                   % 1 THz = 1e12 Hz
dt = 0.05*ps;                 % 时间步长，s
tEnd = 50*ps;                 % 完整时间终点，s
t = (0:round(tEnd/dt))'*dt;   % 时间列，s

t0 = 5*ps;                    % 参考脉冲中心，s
sigma = 4*ps;                % 高斯包络宽度，s
A = 0.6;                      % 样品幅值 / 参考幅值，无量纲
delay = 2*ps;                 % 样品延迟，s
E0 = 1;                       % 模拟电场标度，V/m（人为设定，非实测校准）
f1 = 0.3*THz;                 % 正弦频率，Hz
f2 = 0.6*THz;
f3 = 0.9*THz;

x_start = 0;                  % 时域绘图起始时间，s
x_end = tEnd;                 % 时域绘图结束时间，s

fftStart = 0;                 % FFT 起始时间，s
fftEnd = 20*ps;               % FFT 结束时间，s

%% 2. 生成参考和样品信号
pulse = @(x) ...
    E0 * (sin(2*pi*f1*x) ...
   + 0.5*sin(2*pi*f2*x)  ...
    + 0.2*sin(2*pi*f3*x)) ...
   .* exp(-x.^2/(2*sigma^2));


E_ref = pulse(t - t0);
E_sam = A * pulse(t - t0 - delay);

% % 两条信号使用相同的归一化系数
% scale = max(abs(E_ref));
% E_ref = E_ref / scale;
% E_sam = E_sam / scale;

%% 3. 选择参与 FFT 的部分
timeTolerance = dt*1e-6;     % 避免浮点舍入漏掉窗口端点
idx = t >= fftStart-timeTolerance & t <= fftEnd+timeTolerance;

t_cut = t(idx);
ref_cut = E_ref(idx);
sam_cut = E_sam(idx);

%% 4. 时域图：绿色原始范围，红色 FFT 范围
signals = {E_ref, E_sam};
names = {'Reference', 'Sample'};

figure('Color', 'w');

for k = 1:2
    subplot(2,1,k);
    plot(t, signals{k}, 'b', 'LineWidth', 1);
    hold on;

    % 绿色框：完整原始信号范围
    hGreen = rectangle('Position', ...
        [t(1), -1.5, t(end)-t(1), 3], ...
        'EdgeColor', [0 0.6 0], 'LineWidth', 1.5);

    % 红色框：参与 FFT 的信号范围
    hRed = rectangle('Position', ...
        [fftStart, -1.4, fftEnd-fftStart, 2.8], ...
        'EdgeColor', 'r', 'LineWidth', 1.5);

    xlim([x_start, x_end]);
    ylim([-2 2]);
    xlabel('Time (s)');
    ylabel('Electric field (V/m)');
    title([names{k}, ' : green = original, red = FFT interval']);
    grid on;
end

%% 5. 参考和样品分别做 FFT
N = length(t_cut);

F_ref = fft(ref_cut);
F_sam = fft(sam_cut);

% 单边频率轴：dt 为 s，因此 f 为 Hz
nPositive = floor(N/2) + 1;
f = (0:nPositive-1)' / (N*dt);

% 单边幅度谱
amp_ref = abs(F_ref(1:nPositive)) / N;
amp_sam = abs(F_sam(1:nPositive)) / N;

if mod(N,2) == 0
    amp_ref(2:end-1) = 2*amp_ref(2:end-1);
    amp_sam(2:end-1) = 2*amp_sam(2:end-1);
else
    amp_ref(2:end) = 2*amp_ref(2:end);
    amp_sam(2:end) = 2*amp_sam(2:end);
end

%% 6. 频域图
figure('Color', 'w');
plot(f, amp_ref, 'b', 'LineWidth', 1.2);
hold on;
plot(f, amp_sam, 'r', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('FFT electric field amplitude (V/m)');
legend('Reference', 'Sample');
xlim([0 1*THz]);
grid on;

%% 7. 保存完整原始数据，供你的处理脚本读取
% 三列分别为时间(s)、参考电场(V/m)、样品电场(V/m)
outputFile = 'THz_simulated.xlsx';
outputData = table(t, E_ref, E_sam, 'VariableNames', ...
    {'Time_s', 'E_ref_V_per_m', 'E_sam_V_per_m'});
writetable(outputData, outputFile, 'Sheet', 1, ...
    'WriteMode', 'overwritesheet');
fprintf('已保存：%s（时间 s；参考和样品电场 V/m）\n', outputFile);

