clc; clear;

%% ========================================================================
%  三开口谐振环（三开口环槽）400 组受约束抽样
%
%  模板：三开口谐振环.cst（玻璃衬底 + PEC 环形槽 + 3 条金属桥）
%  依据：课题组的《三开口环槽400组扫参参数设置》文档；下文简称 docx，
%        该文档未随代码发布，相关判据已在脚本里逐条写明
%
%  扫描参数（6 维，单位 um / 度）：p, r, w, g, alpha, beta
%
%  模板几何（读自 .cst 内 Model/3D/Model.mod，已核对）：
%    · 玻璃衬底      solid1，p x p，厚度 t=500，eps_r=4
%    · 环形槽（PEC） solid2，OuterRadius = r，InnerRadius = r-w，厚度 t1=0.2
%    · 金属桥 3 条   solid3/solid4/solid5，宽 g，方位角分别为
%                      0°、90°-alpha、90°+beta
%    · 最终 Solid.Subtract：环形槽 - 3 条桥 = 三段弧槽
%
%  因此三段弧的真实圆心角（本脚本以此为准）：
%      arc1 = 90 + beta
%      arc2 = 360 - alpha - beta      <-- docx 中写为 gamma
%      arc3 = alpha - 90
%  docx 写 "alpha=beta=145 时第三夹角 70°"：360-145-145 = 70，与 arc2 一致，
%  说明 docx 的 gamma 就是 arc2。但 alpha=beta=145 时 arc3 = 55°，才是真正的
%  最短弧；alpha=120 时 arc3 只有 30°。故本脚本默认用真实最短弧做长度约束。
%
%  约束：
%    (1) p - 2*r        >= 40       环外缘到单元边界 >= 20
%    (2) r - w          >= 60       环内侧金属区宽度
%    (3) gamma          >= 55       gamma = 360-alpha-beta
%    (4) L_min          >= 35       L = (r-w/2)*theta*pi/180 - g
%                                   theta 取 STRICT_MIN_ARC 指定的弧
%
%  分组（合计 400，docx 表 4）：
%    主区   p=250~275 : 对称 64  + 非对称 256 = 320
%    边界区 p=275~300 : 对称 16  + 非对称  64 =  80
%  参考点 (260,110,40,40,145,145) 必出现在主区-对称组。
%% ========================================================================

%% ---------------- 开关 ----------------
STRICT_MIN_ARC = true;   % true : 用真实最短弧（arc3/arc2 取小）判 L>=35
                         % false: 完全按 docx 原文，用 arc2(=gamma) 判 L>=35
N_GEN_MULT     = 8;      % 每批生成量 = 需求数 * N_GEN_MULT（冗余倍数）
ROUND_DIGITS   = 1;      % 长度与角度均保留 1 位小数
MIN_ABS_ANGLE  = 1.0;    % 非对称组要求 |alpha-beta| >= 1.0 度（避免舍入后退化成对称）
OUT_DIR = fileparts(mfilename('fullpath'));   % 输出目录：默认与脚本同目录

%% ---------------- 参数范围（docx 表 1） ----------------
p_range     = [250, 300];
r_range     = [100, 130];
w_range     = [ 35,  50];
g_range     = [ 35,  50];
ang_range   = [135, 155];

%% ---------------- 分组（docx 表 4） ----------------
cells = struct( ...
    'name',   {'主区-对称',   '主区-非对称', '边界-对称',   '边界-非对称'}, ...
    'p_lo',   {p_range(1),    p_range(1),  275,           275        }, ...
    'p_hi',   {275,           275,          p_range(2),    p_range(2)}, ...
    'n',      {64,            256,          16,            64         }, ...
    'sym',    {true,          false,        true,          false      });

REF_POINT = [260, 110, 40, 40, 145, 145];   % p r w g alpha beta

rng('shuffle');

%% ---------------- 逐组抽样 ----------------
T_all = zeros(0, 7);          % [groupIdx, p r w g alpha beta]
for c = 1:numel(cells)
    need  = cells(c).n;
    isRef = (c == 1);          % 参考点放主区-对称组
    if isRef, need = need - 1; end

    got = zeros(0, 6);
    tries = 0;
    while size(got, 1) < need && tries < 400
        tries = tries + 1;
        n_gen = max(400, need * N_GEN_MULT);
        batch = sample_batch(n_gen, cells(c).p_lo, cells(c).p_hi, ...
                             cells(c).sym, STRICT_MIN_ARC, MIN_ABS_ANGLE, ...
                             p_range, r_range, w_range, g_range, ang_range);

        if ~isempty(batch)
            batch = round(batch, ROUND_DIGITS);
            batch = batch(validate_all(batch, STRICT_MIN_ARC, cells(c).sym, ...
                             MIN_ABS_ANGLE, p_range, r_range, w_range, g_range, ang_range), :);
            got = [got; batch];                                   %#ok<AGROW>
            if size(got, 1) >= need
                break;
            end
        end
    end

    if size(got, 1) < need
        error('分组 %s 抽样失败：只得到 %d/%d 组。请放宽约束或增大 N_GEN_MULT。', ...
              cells(c).name, size(got, 1), need);
    end

    got = unique(got, 'rows', 'stable');
    got = got(1:need, :);

    if isRef
        got = [REF_POINT; got];
    end

    T_all = [T_all; repmat(c, size(got, 1), 1), got];             %#ok<AGROW>
end

%% ---------------- 汇总 ----------------
ID = arrayfun(@(k) sprintf('run%04d', k), (1:size(T_all,1))', 'UniformOutput', false);
G  = T_all(:, 1);
D  = T_all(:, 2:7);
p = D(:,1); r = D(:,2); w = D(:,3); g = D(:,4); alpha = D(:,5); beta = D(:,6);

gamma = 360 - alpha - beta;
arc1  = 90 + beta;
arc2  = 360 - alpha - beta;
arc3  = alpha - 90;
arcmin = min([arc1, arc2, arc3], [], 2);

L1 = (r - w/2) .* arc1 .* pi/180 - g;
L2 = (r - w/2) .* arc2 .* pi/180 - g;
L3 = (r - w/2) .* arc3 .* pi/180 - g;
Lmin = min([L1, L2, L3], [], 2);

%% ---------------- 写 1：仿真输入表（制表符分隔，6 位小数统一格式） ----------------
fname = fullfile(OUT_DIR, 'SanKaiKou_LHS_400.txt');
fid = fopen(fname, 'w');
fprintf(fid, 'ID\tp\tr\tw\tg\talpha\tbeta\n');
for i = 1:size(D,1)
    fprintf(fid, '%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
        ID{i}, p(i), r(i), w(i), g(i), alpha(i), beta(i));
end
fclose(fid);

%% ---------------- 写 2：审计表（含派生量，便于检查与出图） ----------------
fid = fopen(fullfile(OUT_DIR, 'SanKaiKou_LHS_400_audit.csv'), 'w');
fprintf(fid, 'ID,group,region,class,p,r,w,g,alpha,beta,gamma,arc1,arc2,arc3,arc_min,L1,L2,L3,Lmin,p_minus_2r,r_minus_w\n');
for i = 1:size(D,1)
    fprintf(fid, '%s,%d,%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', ...
        ID{i}, G(i), region_name(G(i)), class_name(alpha(i), beta(i)), ...
        p(i), r(i), w(i), g(i), alpha(i), beta(i), ...
        gamma(i), arc1(i), arc2(i), arc3(i), arcmin(i), ...
        L1(i), L2(i), L3(i), Lmin(i), p(i)-2*r(i), r(i)-w(i));
end
fclose(fid);

%% ---------------- 写 3：试算 32 组（docx：正式 400 组前先试算 24~36 组） ----------------
N_PILOT = 32;
pilotIdx = pick_pilot(D, G, N_PILOT);
fid = fopen(fullfile(OUT_DIR, 'SanKaiKou_pilot32.txt'), 'w');
fprintf(fid, 'ID\tp\tr\tw\tg\talpha\tbeta\n');
for k = 1:numel(pilotIdx)
    i = pilotIdx(k);
    fprintf(fid, '%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
        ID{i}, p(i), r(i), w(i), g(i), alpha(i), beta(i));
end
fclose(fid);

%% ---------------- 控制台报告 ----------------
fprintf('\n==== 已生成 %d 组 ====\n', size(D,1));
fprintf('  分组计数：主区-对称 %d，主区-非对称 %d，边界-对称 %d，边界-非对称 %d\n', ...
    sum(G==1), sum(G==2), sum(G==3), sum(G==4));
fprintf('  参考点 (260,110,40,40,145,145) 已包含：%s\n', ...
    string(any(all(abs(D - REF_POINT) < 1e-9, 2))));
fprintf('  alpha==beta 组数 = %d（docx 要求 80）\n', sum(abs(alpha-beta) < 1e-9));
fprintf('\n  参数          最小      最大      均值\n');
nm = {'p    ','r    ','w    ','g    ','alpha','beta '};
for j = 1:6
    fprintf('  %s  %9.1f %9.1f %9.1f\n', nm{j}, min(D(:,j)), max(D(:,j)), mean(D(:,j)));
end
fprintf('\n  派生量        最小      最大\n');
fprintf('  gamma       %9.1f %9.1f\n', min(gamma), max(gamma));
fprintf('  arc_min     %9.1f %9.1f\n', min(arcmin), max(arcmin));
fprintf('  L_min       %9.1f %9.1f\n', min(Lmin), max(Lmin));
fprintf('  p-2r        %9.1f %9.1f\n', min(p-2*r), max(p-2*r));
fprintf('  r-w         %9.1f %9.1f\n', min(r-w), max(r-w));

% 关键体检：若只按 docx 的 arc2 判据，有多少组其实短弧过短
pass2 = L2 >= 35 - 1e-9;
passmin = Lmin >= 35 - 1e-9;
fprintf('\n  docx(arc2) 判据通过 %d 组；真实最短弧判据通过 %d 组\n', ...
    sum(pass2), sum(passmin));
if STRICT_MIN_ARC
    fprintf('  -> STRICT_MIN_ARC = true，本表已按真实最短弧过滤\n');
else
    fprintf('  -> STRICT_MIN_ARC = false，本表按 docx 原文过滤；其中 %d 组真实最短弧 < 35 um\n', ...
        sum(pass2 & ~passmin));
end
fprintf('\n  已写 SanKaiKou_LHS_400.txt / SanKaiKou_LHS_400_audit.csv / SanKaiKou_pilot32.txt\n');

%% ========================================================================
%                              局部函数
%% ========================================================================
function B = sample_batch(n, p_lo, p_hi, symFlag, strictArc, minDiff, ...
                          p_range, r_range, w_range, g_range, ang_range)
% 顺序条件分层抽样：w -> r -> p -> alpha -> beta -> g
U = lhs_uniform(n, 6);

% 1) w
w = w_range(1) + U(:,1) * (w_range(2) - w_range(1));

% 2) r：下界受 r-w>=60；上界 130
r_lo = w + 60;
r_hi = r_range(2);
ok_r = r_lo <= r_hi;
r = r_lo + U(:,2) .* (r_hi - r_lo);

% 3) p：下界受 p-2r>=40 与本组区间；上界本组区间
p_lo2 = max(2*r + 40, max(p_lo, p_range(1)));
p_hi2 = min(p_hi, p_range(2));
ok_p = p_lo2 <= p_hi2;
p = p_lo2 + U(:,3) .* (p_hi2 - p_lo2);

% 4) alpha / beta
if symFlag
    a_hi = min(ang_range(2), (360 - 55)/2);                % gamma=360-2a >= 55 => a <= 152.5
    alpha = ang_range(1) + U(:,4) * (a_hi - ang_range(1));
    beta  = alpha;
else
    alpha = ang_range(1) + U(:,4) * (ang_range(2) - ang_range(1));
    beta  = ang_range(1) + U(:,5) * (ang_range(2) - ang_range(1));
end
ok_ang = (360 - alpha - beta) >= 55;
if ~symFlag
    ok_ang = ok_ang & (abs(alpha - beta) >= minDiff);
end

% 5) g：由 最短弧净槽长 >= 35 反解上界
arc1 = 90 + beta;
arc2 = 360 - alpha - beta;
arc3 = alpha - 90;
if strictArc
    th = min([arc1, arc2, arc3], [], 2);
else
    th = arc2;
end
g_hi = (r - w/2) .* th .* pi/180 - 35;
g = g_range(1) + U(:,6) .* (g_range(2) - g_range(1));
ok_g = (g <= g_hi);            % 拒绝法：g 保持均匀，避免与 alpha 产生相关性

valid = ok_r & ok_p & ok_ang & ok_g;
B = [p(valid), r(valid), w(valid), g(valid), alpha(valid), beta(valid)];
end

function U = lhs_uniform(n, d)
U = zeros(n, d);
for j = 1:d
    t = ((0:n-1)' + rand(n,1)) / n;
    U(:,j) = t(randperm(n));
end
end

function ok = validate_all(B, strictArc, symFlag, minDiff, ...
                           p_range, r_range, w_range, g_range, ang_range)
p = B(:,1); r = B(:,2); w = B(:,3); g = B(:,4); alpha = B(:,5); beta = B(:,6);
tol = 1e-9;

arc1 = 90 + beta;
arc2 = 360 - alpha - beta;
arc3 = alpha - 90;
if strictArc
    th = min([arc1, arc2, arc3], [], 2);
else
    th = arc2;
end
L = (r - w/2) .* th .* pi/180 - g;

ok = (p >= p_range(1)-tol) & (p <= p_range(2)+tol) & ...
     (r >= r_range(1)-tol) & (r <= r_range(2)+tol) & ...
     (w >= w_range(1)-tol) & (w <= w_range(2)+tol) & ...
     (g >= g_range(1)-tol) & (g <= g_range(2)+tol) & ...
     (alpha >= ang_range(1)-tol) & (alpha <= ang_range(2)+tol) & ...
     (beta  >= ang_range(1)-tol) & (beta  <= ang_range(2)+tol) & ...
     (p - 2*r >= 40 - tol) & ...
     (r - w   >= 60 - tol) & ...
     (360 - alpha - beta >= 55 - tol) & ...
     (L >= 35 - tol);
if symFlag
    ok = ok & (abs(alpha - beta) < 1e-9);
else
    ok = ok & (abs(alpha - beta) >= minDiff - tol);
end
end

function s = region_name(gidx)
if gidx <= 2, s = 'main'; else, s = 'edge'; end
end

function s = class_name(a, b)
if abs(a - b) < 1e-9, s = 'sym'; else, s = 'asym'; end
end

function idx = pick_pilot(D, G, n)
% 试算点选择：先钉住关键锚点，再按分组做 maximin 铺开
N = size(D,1);
groups = unique(G(:));
Z = (D - min(D,[],1)) ./ max(eps, (max(D,[],1) - min(D,[],1)));

refRow = find(all(abs(D - [260 110 40 40 145 145]) < 1e-9, 2), 1);
if isempty(refRow), refRow = 1; end
[~, iRmin] = min(D(:,2));  [~, iRmax] = max(D(:,2));
[~, iPmin] = min(D(:,1));  [~, iPmax] = max(D(:,1));
[~, iAmin] = min(D(:,5));  [~, iAmax] = max(D(:,5));
sel = unique([refRow, iRmin, iRmax, iPmin, iPmax, iAmin, iAmax]', 'stable');
sel = sel(:);

per = ones(numel(groups),1) * floor(n/numel(groups));
per(1:mod(n,numel(groups))) = per(1:mod(n,numel(groups))) + 1;

for gi = 1:numel(groups)
    g = groups(gi);
    pool = find(G == g);
    have = sel(G(sel) == g);
    need = per(gi) - numel(have);
    if need <= 0, continue; end
    Zp = Z(pool,:);
    if isempty(have)
        c = median(Zp,1);
        [~, j0] = min(sum((Zp - c).^2, 2));
        chosen = j0;
        dmin = sqrt(sum((Zp - Zp(j0,:)).^2, 2));
    else
        chosen = zeros(0,1);
        dmin = min(sqrt(sum((Z(pool,:) - Z(have(1),:)).^2, 2)), ...
                   [], 2) * 0 + inf;
        for hh = have(:)'
            dmin = min(dmin, sqrt(sum((Zp - Z(hh,:)).^2, 2)));
        end
    end
    while numel(chosen) < need && numel(chosen) < numel(pool)
        d2 = dmin; d2(chosen) = -inf;
        [~, j] = max(d2);
        chosen(end+1) = j;                                    %#ok<AGROW>
        dmin = min(dmin, sqrt(sum((Zp - Zp(j,:)).^2, 2)));
    end
    add = pool(chosen);
    add = add(~ismember(add, sel));
    sel = [sel(:); add(:)];                                   %#ok<AGROW>
end

sel = unique(sel, 'stable');
idx = sel(1:min(n, numel(sel)));
end
