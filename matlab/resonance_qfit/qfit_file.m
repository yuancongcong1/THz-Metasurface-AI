function result=qfit_file(filename,varargin)
%QFIT_FILE File-to-result convenience entry point with a readable summary.
% r=qfit_file('E:/YCC/20261026/p260r110.txt','InputKind','amplitude', ...
%   'Model','fano','Window',[.466 .471]);
ip=inputParser;
addParameter(ip,'Model','fano');addParameter(ip,'InputKind','amplitude');
addParameter(ip,'Window',[]);addParameter(ip,'Feature','peak');
addParameter(ip,'Baseline','linear');addParameter(ip,'FrequencyUnit','THz');
addParameter(ip,'FrequencyColumn',1);addParameter(ip,'ValueColumn',2);
addParameter(ip,'ImagColumn',[]);addParameter(ip,'NumStarts',12);
addParameter(ip,'OutputPrefix','');addParameter(ip,'ShowFigure',true);
addParameter(ip,'MakePlot',true);parse(ip,varargin{:});cfg=ip.Results;
assert(isfile(filename),'Cannot find spectrum: %s',filename);
assert(xor(strcmpi(cfg.InputKind,'complex'),isempty(cfg.ImagColumn)), ...
    'Set InputKind=complex and ImagColumn together.');
[f,y]=read_cst_spectrum(filename,cfg.FrequencyColumn,cfg.ValueColumn,cfg.ImagColumn);
fprintf('\nData: %s\nFrequency: %.9g to %.9g %s, %d exported points\n', ...
    filename,min(f),max(f),cfg.FrequencyUnit,numel(f));
result=fit_resonance_q(f,y,'Model',cfg.Model,'InputKind',cfg.InputKind, ...
    'Window',cfg.Window,'Feature',cfg.Feature,'Baseline',cfg.Baseline, ...
    'FrequencyUnit',cfg.FrequencyUnit,'NumStarts',cfg.NumStarts);
if isempty(cfg.OutputPrefix)
    [folder,stem]=fileparts(filename);
    prefix=fullfile(folder,'qfit_results',[stem,'_',datestr(now,'yyyymmdd_HHMMSSFFF')]);
else,prefix=char(cfg.OutputPrefix);end
inputCanonical=java.io.File(filename).getCanonicalPath();
for ext={'.csv','.json','.png'}
    targetCanonical=java.io.File([prefix,ext{1}]).getCanonicalPath();
    assert(~strcmpi(char(inputCanonical),char(targetCanonical)),'Output would overwrite input.');
end
save_qfit_result(result,prefix,cfg.MakePlot);
fprintf('\nModel: %s\nPole center f0: %.9g %s\nVisible feature: %.9g %s\n', ...
    result.model,result.f0,result.frequency_unit,result.feature_frequency,result.frequency_unit);
fprintf('Pole Q (Lorentz/Fano): %s\n',format_value(result.q_resonance));
fprintf('Half-contrast spectral Q: %s\n',format_value(result.q_spectral));
fprintf('Absolute half-power / 3 dB Q: %s\n',format_value(result.q_3db));
fprintf('Power-domain RMSE: %.6g\nOutputs: %s.*\n',result.diagnostics.rmse_power,prefix);
for k=1:numel(result.warnings),fprintf('CHECK: %s\n',result.warnings{k});end
fprintf('Q from exported spectra does not verify CST mesh/frequency convergence.\n');
if cfg.ShowFigure
    fig=figure('Color','w','Name','Q fit: inspect spectrum and residuals','NumberTitle','off');
    tiledlayout(fig,2,1);
    nexttile;plot(result.frequency,result.power,'.','DisplayName','Input power');hold on;
    plot(result.frequency,result.fitted_power,'-','DisplayName',result.model);
    plot(result.frequency,result.fitted_baseline,'--','DisplayName','Local baseline');
    ylabel('Power');legend('Location','best');title('Fit: inspect both sides of the resonance');
    nexttile;plot(result.frequency,result.power-result.fitted_power,'.');yline(0,':');
    ylabel('Residual');xlabel(['Frequency / ',result.frequency_unit]);
end
end

function text=format_value(value)
if isempty(value),text='unavailable';else,text=sprintf('%.6g',value);end
end
