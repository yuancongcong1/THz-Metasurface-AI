function summary=batch_qfit(folder,pattern,outputDir,varargin)
%BATCH_QFIT One spectrum per file, one selected mode per fit.
% batch_qfit('D:/spectra','*.txt','D:/fits','Model','fano', ...
%            'InputKind','amplitude','Window',[.30 .36],'Feature','peak')
% f/y columns default to 1/2. For other formats call read_cst_spectrum and
% fit_resonance_q directly. Complex files need explicit column selection.
assert(isfolder(folder),'Input folder not found.');
assert(~strcmp(char(java.io.File(folder).getCanonicalPath()), ...
    char(java.io.File(outputDir).getCanonicalPath())),'Use a separate output directory.');
files=dir(fullfile(folder,pattern));files=files(~[files.isdir]);
assert(~isempty(files),'No matching spectrum files.');
if ~isfolder(outputDir),mkdir(outputDir);end
rows=cell(numel(files),9);
for k=1:numel(files)
    name=files(k).name;rows(k,:)={name,'failed','','','','','','',''};
    try
        [f,y]=read_cst_spectrum(fullfile(folder,name));
        r=fit_resonance_q(f,y,varargin{:});
        save_qfit_result(r,fullfile(outputDir,[name,'_fit']));
        rows(k,:)={name,'fit_completed',r.model,num2str(r.f0,12),num2str(r.q_resonance,12), ...
            num2str(r.q_spectral,12),num2str(r.diagnostics.rmse_power,12),strjoin(r.warnings,'; '),''};
    catch err
        rows{k,9}=err.message;
    end
    fprintf('%s %s\n',rows{k,2},name);
end
summary=table(rows(:,1),rows(:,2),rows(:,3),rows(:,4),rows(:,5),rows(:,6),rows(:,7),rows(:,8),rows(:,9), ...
    'VariableNames',{'file','status','model','f0','q_resonance','q_spectral','rmse_power','warnings','error'});
writetable(summary,fullfile(outputDir,'summary.csv'));
end
