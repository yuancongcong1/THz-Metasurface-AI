%% Run a parameter table through a CST template via Windows COM automation.
% Edit CONFIGURATION for your model. Every numeric table column except ID
% must match a CST parameter name. MATLAB + CST must run on the same Windows PC.

%% CONFIGURATION
cfg.inputFile = fullfile(pwd, 'parameters.tsv');
cfg.templateFile = fullfile(pwd, 'model.cst');
cfg.workDir = fullfile(pwd, 'cst_runs');
cfg.startIndex = 1;
cfg.stopIndex = inf;
cfg.skipExisting = true;
cfg.continueOnError = true;
cfg.resultWaitSec = 120;   % Wait for exports AFTER solver.Start returns.
cfg.pollSec = 2;
cfg.solverObject = 'FDSolver';
cfg.solverMethod = 'Start';

% Check these EXACT names in your CST navigation tree. Choose the format
% that the selected CST trace actually exports: 'db20' converts dB to
% linear magnitude; 'linear' keeps the values; 'raw' keeps the two columns.
cfg.results = struct( ...
    'name', {'S21','S11'}, ...
    'tree', {'1D Results\S-Parameters\SZmin(2),Zmax(2)', ...
             '1D Results\S-Parameters\SZmax(2),Zmax(2)'}, ...
    'format', {'db20','db20'});

runSweep(cfg);

%% LOCAL FUNCTIONS
function runSweep(cfg)
%% INPUT CHECKS
assert(isfile(cfg.inputFile), 'Parameter file not found: %s', cfg.inputFile);
assert(isfile(cfg.templateFile), 'CST template not found: %s', cfg.templateFile);
assert(cfg.resultWaitSec > 0 && cfg.pollSec > 0);
assert(~isempty(cfg.results), 'Configure at least one result trace.');
T = readtable(cfg.inputFile, 'FileType', 'text', 'Delimiter', '\t', ...
              'VariableNamingRule', 'preserve');
assert(any(strcmp(T.Properties.VariableNames, 'ID')), 'Input must contain ID.');
assert(height(T) > 0, 'Input table is empty.');
paramNames = setdiff(T.Properties.VariableNames, {'ID'}, 'stable');
assert(~isempty(paramNames), 'Input has no parameter columns.');
ids = string(T.ID);
assert(all(~ismissing(ids) & strlength(ids) > 0), 'IDs cannot be empty.');
assert(numel(unique(ids)) == numel(ids), 'IDs must be unique.');
assert(all(~cellfun('isempty', regexp(cellstr(ids), ...
       '^[A-Za-z0-9_-]+$', 'once'))), 'IDs must be file-safe.');
for k = 1:numel(paramNames)
    x = T.(paramNames{k});
    assert(isnumeric(x) && isvector(x) && all(isfinite(x)), ...
           'Parameter %s must contain finite numeric values.', paramNames{k});
end
resultNames = {cfg.results.name};
assert(numel(unique(resultNames)) == numel(resultNames), ...
       'Result names must be unique.');
for k = 1:numel(cfg.results)
    assert(~isempty(regexp(cfg.results(k).name, '^[A-Za-z0-9_-]+$', 'once')));
    assert(ismember(cfg.results(k).format, {'db20','linear','raw'}));
end
lastIndex = min(height(T), cfg.stopIndex);
assert(cfg.startIndex >= 1 && cfg.startIndex <= lastIndex && ...
       cfg.startIndex == floor(cfg.startIndex), 'Invalid start/stop index.');

resultRoot = fullfile(cfg.workDir, 'results');
projectRoot = fullfile(cfg.workDir, 'projects');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end
if ~exist(projectRoot, 'dir'), mkdir(projectRoot); end
failureLog = fullfile(cfg.workDir, 'failed_runs.tsv');

%% CST SESSION
cst = actxserver('CSTStudio.Application');
appCleanup = onCleanup(@() safeDelete(cst)); %#ok<NASGU>
completed = 0;
skipped = 0;
failed = 0;
for row = cfg.startIndex:lastIndex
    runID = char(ids(row));
    resultDir = fullfile(resultRoot, runID);
    marker = fullfile(resultDir, 'complete.txt');
    if cfg.skipExisting && isComplete(marker, resultDir, cfg.results)
        skipped = skipped + 1;
        fprintf('[%d/%d] %s already complete; skipped.\n', row, lastIndex, runID);
        continue;
    end
    if ~exist(resultDir, 'dir'), mkdir(resultDir); end
    if isfile(marker), delete(marker); end

    % A fresh attempt directory prevents old CST results from being reused.
    attemptDir = fullfile(projectRoot, runID, ...
        ['attempt_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
    if ~exist(attemptDir, 'dir'), mkdir(attemptDir); end
    runFile = fullfile(attemptDir, 'run.cst');
    mws = [];
    try
        copyfile(cfg.templateFile, runFile);
        mws = invoke(cst, 'NewMWS');
        invoke(mws, 'OpenFile', runFile);
        for k = 1:numel(paramNames)
            value = T.(paramNames{k})(row);
            invoke(mws, 'StoreParameter', paramNames{k}, sprintf('%.15g', value));
        end
        invoke(mws, 'Rebuild');
        fprintf('[%d/%d] Simulating %s...\n', row, lastIndex, runID);
        solver = invoke(mws, cfg.solverObject);
        invoke(solver, cfg.solverMethod);

        for k = 1:numel(cfg.results)
            spec = cfg.results(k);
            finalFile = fullfile(resultDir, [runID '_' spec.name '.tsv']);
            tempFile = [finalFile '.tmp'];
            exportResult(mws, spec, tempFile, cfg.resultWaitSec, cfg.pollSec);
            movefile(tempFile, finalFile, 'f');
        end
        invoke(mws, 'Save');
        safeClose(mws);
        mws = [];
        fid = fopen(marker, 'w');
        assert(fid >= 0, 'Could not create completion marker.');
        fprintf(fid, 'ID\t%s\nCompleted\t%s\n', runID, char(datetime('now')));
        fclose(fid);
        completed = completed + 1;
        fprintf('Completed %s.\n', runID);
    catch ME
        failed = failed + 1;
        safeClose(mws);
        if isfile(marker), delete(marker); end
        msg = regexprep(ME.message, '[\r\n\t]+', ' ');
        fid = fopen(failureLog, 'a');
        if fid >= 0
            fprintf(fid, '%s\t%s\n', runID, msg);
            fclose(fid);
        end
        warning('Run %s failed: %s', runID, msg);
        if ~cfg.continueOnError, rethrow(ME); end
    end
end
fprintf('Sweep finished: %d completed, %d skipped, %d failed.\n', ...
        completed, skipped, failed);
end

function exportResult(mws, spec, outputFile, timeoutSec, pollSec)
    t0 = tic;
    lastError = '';
    while toc(t0) < timeoutSec
        try
            if isfile(outputFile), delete(outputFile); end
            invoke(mws, 'SelectTreeItem', spec.tree);
            exporter = invoke(mws, 'ASCIIExport');
            invoke(exporter, 'Reset');
            invoke(exporter, 'FileName', outputFile);
            invoke(exporter, 'Execute');
            assert(isfile(outputFile), 'CST did not create the export file.');
            data = readmatrix(outputFile, 'FileType', 'text');
            assert(size(data,2) >= 2, 'Expected two numeric export columns.');
            data = data(:,1:2);
            data = data(all(isfinite(data),2),:);
            assert(size(data,1) >= 2, 'Export has fewer than two data rows.');
            assert(all(diff(data(:,1)) > 0), ...
                   'Frequency values must be strictly increasing.');
            if strcmp(spec.format, 'db20')
                data(:,2) = 10.^(data(:,2)/20);
            end
            writematrix(data, outputFile, 'Delimiter', 'tab');
            return;
        catch ME
            lastError = ME.message;
            pause(pollSec);
        end
    end
    error('Export %s timed out after %.0f s: %s', ...
          spec.name, timeoutSec, lastError);
end

function tf = isComplete(marker, folder, specs)
    tf = isfile(marker);
    for k = 1:numel(specs)
        files = dir(fullfile(folder, ['*_' specs(k).name '.tsv']));
        tf = tf && isscalar(files) && files(1).bytes > 0;
    end
end

function safeClose(mws)
    if ~isempty(mws)
        try
            invoke(mws, 'Quit');
        catch
        end
    end
end

function safeDelete(cst)
    try
        delete(cst);
    catch
    end
end

