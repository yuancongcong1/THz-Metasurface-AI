clc;
clear;
close all;

%% =========================
% MATLAB - CST 三开口谐振环（三开口环槽）LHS 参数扫描
% 模板：同目录下的 三开口谐振环.cst（*.cst 不入库，需自备干净副本）
% 参数：p, r, w, g, alpha, beta （6 维，um / deg）
% 导出：S21（透射）+ S11（反射），dB -> 线性幅度
% =========================
%
% 相对早期版本改掉的地方（均为实测踩坑）：
%   1) 参数表不再写死绝对路径，统一 fullfile(basePath, ...)，试算/正式只改一行
%   2) 等待与导出合成一个 export_sparam，删掉 60 s 二段等待与两份重复导出代码
%   3) 结果先写 <name>.part.txt，校验通过才改名
%      -> 中断、超时、重跑都不会再把之前算好的结果删掉
%   4) 单组出错只跳过该组并把 CST 实例重开一次，400 组不会因一组崩掉整晚
%   5) 跳过条件改为 S21、S11 都在（旧版只看 S21，缺 S11 会被误判为已完成）
%   6) 进度 / ETA、频点密度体检、sweep_log.csv 续扫对账

%% =========================================
%               配置区（改这里）
%% =========================================

% 路径：默认都放在本脚本所在目录下，换目录改这两行即可
basePath = fileparts(mfilename('fullpath'));

template = fullfile(basePath, '三开口谐振环.cst');

runfile = fullfile(basePath, 'Temp', 'run.cst');

resultPath = fullfile(basePath, 'SanKaiKou_Results');

% ---- 参数表：二选一，跑正式 400 组时把上面那行注释掉 ----
% paramFile = fullfile(basePath, 'SanKaiKou_pilot32.txt');   % 试算 32 组
paramFile = fullfile(basePath, 'SanKaiKou_LHS_400.txt'); % 正式 400 组

startIdx = 1;               % 断点续扫：从第几组开始（例如已完成 300 组可填 301）

skipExisting = true;        % true = S21、S11 都已存在的结果自动跳过（续扫建议 true）

waitTimeoutSec = 3600;      % 单组等结果的最长时间（秒），超时跳过该组
                            % 频域 + 四面体自适应网格，0.1~0.5 THz 单组约几十秒，
                            % 若你机器上明显更久，把这里调大

minPoints = 500;            % 频点少于这个数就提示（提 Q 会不准）

% ---- CST 导航树中的两条曲线 ----
% 名字必须与树上完全一致：SelectTreeItem 对不存在的名字不报错、只静默不切换视图，
% 于是 ASCIIExport 会一直导出当前那条曲线（这正是之前 S11 与 S21 内容一模一样的原因）。
S21_TREE = '1D Results\S-Parameters\SZmin(1),Zmax(1)';   % 透射
S11_TREE = '1D Results\S-Parameters\SZmax(1),Zmax(1)';   % 反射

%% =========================================
%            检查输入与建目录
%% =========================================

if ~exist(template, 'file')
    error('找不到 CST 模板：%s', template);
end
if ~exist(paramFile, 'file')
    error('找不到参数文件：%s', paramFile);
end

if ~exist(resultPath, 'dir')
    mkdir(resultPath);
end
if ~exist(fullfile(basePath, 'Temp'), 'dir')
    mkdir(fullfile(basePath, 'Temp'));
end

%% =========================================
%           读取 LHS 参数文件
%% =========================================

paramData = readtable(paramFile, 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
N = height(paramData);

needVars = {'ID', 'p', 'r', 'w', 'g', 'alpha', 'beta'};
missing = needVars(~ismember(needVars, paramData.Properties.VariableNames));
if ~isempty(missing)
    error('参数文件缺少列：%s', strjoin(missing, ', '));
end

if startIdx < 1 || startIdx > N
    error('startIdx 必须介于 1 和 %d 之间。', N);
end

fprintf('参数文件 : %s\n', paramFile);
fprintf('共 %d 组，从第 %d 组开始\n', N, startIdx);

%% =========================================
%              启动 CST
%% =========================================

logFile = fullfile(resultPath, 'sweep_log.csv');
if ~exist(logFile, 'file')
    fid = fopen(logFile, 'w');
    fprintf(fid, 'ID,p,r,w,g,alpha,beta,status,seconds,S21,S11\n');
    fclose(fid);
end

fprintf('正在启动 CST ...\n');
cst = actxserver('CSTStudio.Application');

%% =========================================
%              开始扫描
%% =========================================

tAll  = tic;
nDone = 0;
nSkip = 0;
nFail = 0;

for i = startIdx:N

    runID = strtrim(char(string(paramData.ID(i))));
    p  = paramData.p(i);
    r  = paramData.r(i);
    w  = paramData.w(i);
    g  = paramData.g(i);
    al = paramData.alpha(i);
    be = paramData.beta(i);

    S21_txt = fullfile(resultPath, [runID '_S21.txt']);
    S11_txt = fullfile(resultPath, [runID '_S11.txt']);

    % ---- 跳过已完成（两条曲线都要在）----
    if skipExisting && has_data(S21_txt) && has_data(S11_txt)
        nSkip = nSkip + 1;
        fprintf('[%d/%d] %s 已存在，跳过\n', i, N, runID);
        continue;
    end

    % ---- 进度与耗时估计 ----
    fprintf('\n=========================\n');
    fprintf('[%d/%d] 开始 %s\n', i, N, runID);
    fprintf('  p=%.2f  r=%.2f  w=%.2f  g=%.2f  alpha=%.2f  beta=%.2f\n', ...
            p, r, w, g, al, be);
    if nDone > 0
        perRun = toc(tAll) / nDone;
        fprintf('  已用 %.1f min，平均 %.2f min/组，预计剩余 %.1f min\n', ...
                toc(tAll)/60, perRun/60, perRun*(N - i + 1)/60);
    end
    fprintf('=========================\n');

    tRun   = tic;
    status = 'ok';

    try
        %% ---- 删除旧 run.cst 与缓存 ----
        if exist(runfile, 'file')
            delete(runfile);
        end
        cacheDir = fullfile(basePath, 'Temp', 'run');
        if exist(cacheDir, 'dir')
            rmdir(cacheDir, 's');
        end

        %% ---- 复制模板并打开 ----
        copyfile(template, runfile);
        mws = invoke(cst, 'NewMWS');
        invoke(mws, 'OpenFile', runfile);
        pause(2);

        %% ---- 写入 6 个扫描参数 ----
        invoke(mws, 'StoreParameter', 'p',     num2str(p,  '%.2f'));
        invoke(mws, 'StoreParameter', 'r',     num2str(r,  '%.2f'));
        invoke(mws, 'StoreParameter', 'w',     num2str(w,  '%.2f'));
        invoke(mws, 'StoreParameter', 'g',     num2str(g,  '%.2f'));
        invoke(mws, 'StoreParameter', 'alpha', num2str(al, '%.2f'));
        invoke(mws, 'StoreParameter', 'beta',  num2str(be, '%.2f'));

        invoke(mws, 'Rebuild');
        pause(2);

        %% ---- 开始仿真 ----
        disp('开始仿真...');
        fdsolver = invoke(mws, 'FDSolver');
        invoke(fdsolver, 'Start');

        %% ---- 等待并导出：能读出数值就说明求解结束 ----
        [ok21, nP21] = export_sparam(mws, S21_TREE, S21_txt, waitTimeoutSec);
        if ~ok21
            status = 'timeout';
        else
            disp('仿真结束，S21 导出完成');
            [ok11, nP11] = export_sparam(mws, S11_TREE, S11_txt, 300);
            if ok11
                fprintf('S21 / S11 导出完成（%d / %d 个频点）\n', nP21, nP11);
            else
                status = 'fail S11';
            end

            % ---- 频点密度体检：太稀就提不出 Q ----
            if nP21 < minPoints
                warning(['本组只导出 %d 个频点，过于稀疏，无法可靠提取 Q。' ...
                         '请在模板里把 Frequency samples 设为 Equidistant 且足够密'], nP21);
            end
        end

        %% ---- 保存并关闭 ----
        invoke(mws, 'Save');
        invoke(mws, 'Quit');
        pause(2);

    catch ME
        status = ['error: ' ME.message];
        warning('第 %d 组 %s 失败：%s', i, runID, ME.message);
        % 工程可能停在半开状态：关掉 CST 再重开，避免污染后面所有组
        try
            invoke(mws, 'Quit');
        catch
        end
        try
            delete(cst);
        catch
        end
        try
            cst = actxserver('CSTStudio.Application');
        catch
        end
        pause(5);
    end

    %% ---- 记日志 ----
    els = toc(tRun);
    if strcmp(status, 'ok')
        nDone = nDone + 1;
    else
        nFail = nFail + 1;
    end

    fid = fopen(logFile, 'a');
    fprintf(fid, '%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%s,%.1f,%d,%d\n', ...
        runID, p, r, w, g, al, be, status, els, ...
        exist(S21_txt, 'file') > 0, exist(S11_txt, 'file') > 0);
    fclose(fid);

    fprintf('本轮结束：%s，用时 %.0f s\n', status, els);
end

%% =========================================
%              关闭 CST
%% =========================================

try
    delete(cst);
catch
end
clear cst;

el = toc(tAll);
fprintf('\n===== 三开口谐振环 LHS 扫描结束 =====\n');
fprintf('本次实测成功 : %d 组\n', nDone);
fprintf('跳过已存在   : %d 组\n', nSkip);
fprintf('失败 / 超时  : %d 组\n', nFail);
fprintf('总耗时       : %.1f min\n', el/60);
fprintf('结果目录     : %s\n', resultPath);

%% ========================================================================
%                              局部函数
%% ========================================================================

function tf = has_data(f)
% 结果文件存在且像样（>200 字节，排除空壳与半截文件）
tf = false;
if exist(f, 'file')
    d = dir(f);
    tf = ~isempty(d) && d.bytes > 200;
end
end

function [ok, nPts] = export_sparam(mws, treeName, outTxt, waitSec)
% 从 CST 导航树导出 S 参数 → 两列 [频率(THz), 线性幅度]
%
% 为什么用"能不能导出数值"判断求解是否结束：
%   CST 2025 的 SelectTreeItem 对不存在或还没算完的曲线都不报错，只是静默不切换
%   视图，所以 try/catch 包住它判断不出求解状态。模板副本里本来没有旧结果，
%   因此一旦能读出两列数字，导出的必然是本次求解的曲线。
%
% 为什么先写 .part.txt 再改名：
%   直接写 outTxt 的话，重跑时会在导出成功前就把上次算好的结果删掉，中途中断或
%   超时就白白丢一组数据。校验通过后 movefile 一次性落盘，旧结果只在有新数据时才被替换。
%
% 中间文件必须也是 .txt：readmatrix 按扩展名选解析器，.raw 之类会直接报错。
ok   = false;
nPts = 0;
part = [outTxt '.part.txt'];

try
    ascii = invoke(mws, 'ASCIIExport');
catch
    return;
end

t0    = tic;
nLoop = 0;

while true
    if exist(part, 'file')
        delete(part);
    end

    % ---- 选曲线 → 导出（算没算完都走一遍，失败就下一轮）----
    try
        invoke(mws, 'SelectTreeItem', treeName);
        invoke(ascii, 'Reset');
        invoke(ascii, 'FileName', part);
        invoke(ascii, 'Execute');
    catch
        % 求解仍在进行或视图未就绪
    end

    % ---- 读回校验 ----
    A = [];
    if exist(part, 'file')
        try
            A = readmatrix(part, 'FileType', 'text');
        catch
            % 文件还没写完整
        end
    end

    if ~isempty(A) && size(A, 2) >= 2
        freq = A(:, 1);                 % THz
        mag  = 10 .^ (A(:, 2) / 20);    % CST ASCII 导出第二列是 dB
        writematrix([freq, mag], part, 'Delimiter', 'tab');
        movefile(part, outTxt, 'f');    % 校验通过才覆盖旧结果
        ok   = true;
        nPts = numel(freq);
        if max(mag) > 1.0001
            warning('%s 的 |S| 最大 %.4f > 1，第二列可能不是 dB', outTxt, max(mag));
        end
        return;
    end

    % ---- 超时 ----
    if toc(t0) > waitSec
        warning('等待 %d s 仍导不出 %s（求解未结束或曲线名不对）', waitSec, treeName);
        if exist(part, 'file')
            delete(part);
        end
        return;
    end

    nLoop = nLoop + 1;
    if mod(nLoop, 10) == 0
        fprintf('  等待仿真结果... (%.0f s)\n', toc(t0));
    end
    pause(3);
end
end
