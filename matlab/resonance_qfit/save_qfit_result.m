function save_qfit_result(result,prefix,makePlot)
%SAVE_QFIT_RESULT Export metadata JSON, power fit/residual CSV and PNG.
% Set makePlot=false on machines without MATLAB graphics support.
if nargin<3,makePlot=true;end
parent=fileparts(prefix);if ~isempty(parent) && ~isfolder(parent),mkdir(parent);end
metadata=rmfield(result,{'frequency','power','fitted_power','fitted_baseline'});
fid=fopen([prefix,'.json'],'w','n','UTF-8');assert(fid>=0,'Cannot create JSON file.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(metadata,'PrettyPrint',true));clear cleanup;
data=table(result.frequency,result.power,result.fitted_power,result.fitted_baseline, ...
    result.power-result.fitted_power,'VariableNames',{'frequency','power','fitted_power','baseline','residual'});
writetable(data,[prefix,'.csv']);
if ~makePlot,return;end
fig=figure('Visible','off','Color','w');guard=onCleanup(@() close(fig));
tiledlayout(2,1);
nexttile;
plot(result.frequency,result.power,'.','DisplayName','Input power');hold on;
plot(result.frequency,result.fitted_power,'-','DisplayName',result.model);
plot(result.frequency,result.fitted_baseline,'--','DisplayName','Local baseline');
xline(result.f0,':','f0','HandleVisibility','off');ylabel('Power response');legend('Location','best');
title(sprintf('%s  pole Q=%s  spectral Q=%s',result.model, ...
    mat2str(result.q_resonance,5),mat2str(result.q_spectral,5)),'Interpreter','none');
nexttile;plot(result.frequency,result.power-result.fitted_power,'.');yline(0,':');
xlabel(['Frequency / ',result.frequency_unit]);ylabel('Residual');
exportgraphics(fig,[prefix,'.png'],'Resolution',180);
end
