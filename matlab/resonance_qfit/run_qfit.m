% Select one CST TXT/CSV/TSV export and fit one isolated resonance.
% Change these settings to match your exported data. No toolboxes required.
cfg.model='fano';             % auto, lorentzian, fano, gaussian, voigt, pseudo_voigt
cfg.inputKind='amplitude';    % amplitude, power, db, complex
cfg.feature='peak';           % peak, dip, auto
cfg.baseline='linear';        % constant or linear
cfg.frequencyUnit='THz';
cfg.fitWindow=[];             % e.g. [0.30 0.35]; empty opens a range dialog
cfg.fcol=1;cfg.ycol=2;cfg.imagcol=[]; % MATLAB column numbers are ONE-based

[name,folder]=uigetfile({'*.txt;*.csv;*.tsv;*.dat','Spectrum TXT/CSV/TSV/DAT'});
if isequal(name,0),return;end
if xor(strcmp(cfg.inputKind,'complex'),~isempty(cfg.imagcol))
    error('Set InputKind=complex and imagcol together.');
end
[f,y]=read_cst_spectrum(fullfile(folder,name),cfg.fcol,cfg.ycol,cfg.imagcol);
if isempty(cfg.fitWindow)
    answer=inputdlg({'Lower fit frequency','Upper fit frequency'}, ...
        'Isolate one resonance (same units as the file)',[1 45],{num2str(min(f)),num2str(max(f))});
    if isempty(answer),return;end
    cfg.fitWindow=[str2double(answer{1}),str2double(answer{2})];
end
result=fit_resonance_q(f,y,'Model',cfg.model,'InputKind',cfg.inputKind, ...
    'Feature',cfg.feature,'Baseline',cfg.baseline,'Window',cfg.fitWindow, ...
    'FrequencyUnit',cfg.frequencyUnit);
disp(rmfield(result,{'frequency','power','fitted_power','fitted_baseline'}));
for k=1:numel(result.warnings),warning('%s',result.warnings{k});end
[~,stem]=fileparts(name);
prefix=fullfile(folder,'qfit_results',[stem,'_',datestr(now,'yyyymmdd_HHMMSS')]);
save_qfit_result(result,prefix);
fprintf('Saved fit, residuals and metadata: %s\n',prefix);
