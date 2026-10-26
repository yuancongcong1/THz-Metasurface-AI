function START_HERE
% A task menu for this repository. Open in MATLAB and press Run.
repoRoot=fileparts(mfilename('fullpath'));
choice=menu('超表面代码仓库：请选择要做的事', ...
    '计算 Q 值：选择频谱和拟合窗口', ...
    '计算 Q 值：运行 p260r110 数据示例', ...
    '计算 Q 值：对比 p260r110 两个偏振', ...
    '生成一份模拟 THz 时域数据', ...
    '打开实际 THz 时域数据处理设置', ...
    '打开 CST 参数采样设置', ...
    '打开 MATLAB-CST 批量仿真设置', ...
    '打开第一次使用指南');
switch choice
    case 1
        addpath(fullfile(repoRoot,'matlab','resonance_qfit'));qfit_wizard;
    case 2
        run(fullfile(repoRoot,'EXAMPLE_P260R110.m'));
    case 3
        run(fullfile(repoRoot,'EXAMPLE_P260R110_TWO_POL.m'));
    case 4
        addpath(fullfile(repoRoot,'matlab','tds'));
        outputFolder=fullfile(repoRoot,'results','tds_demo');
        if ~isfolder(outputFolder),mkdir(outputFolder);end
        fprintf('Demo output folder: %s\n',outputFolder);
        previousFolder=pwd;cd(outputFolder);
        % The script clears its own workspace; preserve cleanup outside it.
        try,evalin('base','data_sin_shengcheng');catch err,cd(previousFolder);rethrow(err);end
        cd(previousFolder);
    case 5
        edit(fullfile(repoRoot,'matlab','tds','shujuchuli.m'));
        fprintf('Set time/reference/sample columns, units and FFT interval first. See docs/START_HERE.md.\n');
    case 6
        edit(fullfile(repoRoot,'matlab','cst','generate_parameters.m'));
        fprintf('This is a geometry sampling template. Set the parameters and constraints for YOUR structure first.\n');
    case 7
        edit(fullfile(repoRoot,'matlab','cst','run_cst_sweep.m'));
        fprintf('Configure CST project, parameter names and result-tree path before running.\n');
    case 8
        edit(fullfile(repoRoot,'docs','START_HERE.md'));
end
end
