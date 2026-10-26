% Run both exported polarization channels. Their x/y naming is not needed.
% Open this file in MATLAB and press Run. Select the data folder if needed.
repoRoot=fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'matlab','resonance_qfit'));
folder='E:\YCC\20261026';
if ~isfile(fullfile(folder,'p260r110.txt')) || ...
   ~isfile(fullfile(folder,'p260r110-2.txt'))
    folder=uigetdir(pwd,'Select the folder containing both p260r110 TXT files');
    if isequal(folder,0),return;end
end
fileA=fullfile(folder,'p260r110.txt');
fileB=fullfile(folder,'p260r110-2.txt');
assert(isfile(fileA) && isfile(fileB),'Both p260r110.txt and p260r110-2.txt are required.');
fitWindow=[.468 .482];
fprintf('\n===== Polarization A (mode 1) =====\n');
resultA=qfit_file(fileA,'InputKind','amplitude','Model','fano', ...
    'Feature','peak','Window',fitWindow,'FrequencyUnit','THz');
fprintf('\n===== Polarization B (mode 2) =====\n');
resultB=qfit_file(fileB,'InputKind','amplitude','Model','fano', ...
    'Feature','peak','Window',fitWindow,'FrequencyUnit','THz');
fprintf('\nA: f0=%.9g THz, Q_pole=%.6g, Q_spectral=%.6g\n', ...
    resultA.f0,resultA.q_resonance,resultA.q_spectral);
fprintf('B: f0=%.9g THz, Q_pole=%.6g, Q_spectral=%.6g\n', ...
    resultB.f0,resultB.q_resonance,resultB.q_spectral);
fprintf('The fit describes exported data; verify CST solve points, mesh and glass loss.\n');
