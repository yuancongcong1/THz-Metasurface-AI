% Real-data walkthrough. The input data itself is NOT included in the repo.
% Open this file in MATLAB and press Run. Select p260r110.txt if E: is absent.
repoRoot=fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'matlab','resonance_qfit'));
inputFile='E:\YCC\20261026\p260r110.txt';
if ~isfile(inputFile)
    [name,folder]=uigetfile({'*.txt','CST magnitude TXT'},'Select your p260r110.txt');
    if isequal(name,0),return;end
    inputFile=fullfile(folder,name);
end
result=qfit_file(inputFile,'InputKind','amplitude','Model','fano', ...
    'Feature','peak','Window',[.466 .471],'FrequencyUnit','THz');
fprintf('\nThis example window applies only to the supplied p260r110 spectrum.\n');
fprintf('For another spectrum run START_QFIT and select its own resonance window.\n');
