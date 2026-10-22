%% Generate constrained, space-filling parameter samples for a CST model.
% Edit the CONFIGURATION section for a new geometry. The AuSRR rules below
% are an example, not universal physical or manufacturing constraints.
% MATLAB R2019b or newer (local functions in scripts); no toolbox needed.

%% CONFIGURATION
cfg.sampleCount = 400;
cfg.poolMultiplier = 3;          % Candidate pool size / requested samples
cfg.batchSize = 800;
cfg.maxBatches = 30;             % Fail safely when constraints are impossible
cfg.seed = 20261022;            % Use 'shuffle' for a different run each time
cfg.decimalPlaces = 1;
cfg.selection = 'maximin';       % 'maximin' for coverage; 'random' for random subset
cfg.outputFile = fullfile(pwd, 'parameters.tsv');

% Sampling order must respect dependencies in the bound functions.
cfg.sampleOrder = {'w','d','r','l2','l1','p'};
cfg.outputOrder = {'p','r','w','d','l1','l2'};

% Each function receives already sampled columns in struct s and row count n.
% Bounds may be scalar or n-by-1; all units here are micrometres.
cfg.lower = {
    @(s,n) 20, ...                                      % w
    @(s,n) 20, ...                                      % d
    @(s,n) max(90, s.d/2 + 4*s.w), ...                  % r
    @(s,n) max(40, 2*s.w), ...                          % l2
    @(s,n) max(40, s.l2), ...                           % l1
    @(s,n) max(200, 2*s.r + 20) ...                    % p
};
cfg.upper = {
    @(s,n) 30, ...                                      % w
    @(s,n) 60, ...                                      % d
    @(s,n) 140, ...                                     % r
    @(s,n) 130, ...                                     % l2
    @(s,n) min(170, 2*sqrt(max(0, (s.r-s.w).^2 ...     % l1
                             - (s.d/2+s.w).^2))), ...
    @(s,n) 300 ...                                      % p
};

% Optional extra vectorized constraint, evaluated after rounding. Return
% one logical value per row. The example is already covered by the bounds.
cfg.extraConstraint = @(s) true(size(s.p));

% Fixed, rough ranges used ONLY to scale distances for maximin selection.
% They do not override the conditional bounds above.
cfg.distanceRanges = [200 300; 90 140; 20 30; 20 60; 40 170; 40 130];

%% GENERATION
assert(cfg.sampleCount >= 1 && cfg.batchSize >= 1 && cfg.maxBatches >= 1);
assert(numel(cfg.sampleOrder) == numel(cfg.lower) && ...
       numel(cfg.sampleOrder) == numel(cfg.upper));
assert(isequal(sort(cfg.sampleOrder), sort(cfg.outputOrder)), ...
       'Output names must be a permutation of sampling names.');
assert(size(cfg.distanceRanges,1) == numel(cfg.outputOrder) && ...
       all(cfg.distanceRanges(:,2) > cfg.distanceRanges(:,1)));
assert(ismember(cfg.selection, {'maximin','random'}));
rng(cfg.seed);

poolTarget = max(cfg.sampleCount, ceil(cfg.sampleCount*cfg.poolMultiplier));
pool = zeros(0, numel(cfg.outputOrder));
generated = 0;
feasible = 0;
for batchNo = 1:cfg.maxBatches
    n = cfg.batchSize;
    generated = generated + n;
    U = lhsUnit(n, numel(cfg.sampleOrder));
    s = struct();
    valid = true(n,1);
    for j = 1:numel(cfg.sampleOrder)
        lo = expandBound(cfg.lower{j}(s,n), n);
        hi = expandBound(cfg.upper{j}(s,n), n);
        valid = valid & isfinite(lo) & isfinite(hi) & hi >= lo;
        s.(cfg.sampleOrder{j}) = lo + U(:,j).*(hi-lo);
    end
    s = roundStruct(s, cfg.sampleOrder, cfg.decimalPlaces);
    valid = valid & checkBounds(s, cfg, n) & ...
            logical(expandBound(cfg.extraConstraint(s), n));
    feasible = feasible + nnz(valid);
    block = zeros(nnz(valid), numel(cfg.outputOrder));
    for j = 1:numel(cfg.outputOrder)
        x = s.(cfg.outputOrder{j});
        block(:,j) = x(valid);
    end
    pool = unique([pool; block], 'rows', 'stable');
    fprintf('Batch %d: %d feasible, %d unique in pool / %d target\n', ...
            batchNo, nnz(valid), size(pool,1), poolTarget);
    if size(pool,1) >= poolTarget
        break;
    end
end
if size(pool,1) < cfg.sampleCount
    error('Only %d unique feasible samples after %d batches. Check bounds.', ...
          size(pool,1), cfg.maxBatches);
end
if size(pool,1) < poolTarget
    warning('Candidate pool is smaller than requested; using %d rows.', size(pool,1));
end

if strcmp(cfg.selection, 'maximin')
    chosen = farthestPointSubset(pool, cfg.sampleCount, cfg.distanceRanges);
else
    chosen = randperm(size(pool,1), cfg.sampleCount);
end
samples = pool(chosen,:);
samples = samples(randperm(cfg.sampleCount),:); % Randomize run order.
ids = compose('run%04d', (1:cfg.sampleCount)');
T = [table(ids, 'VariableNames', {'ID'}), ...
     array2table(samples, 'VariableNames', cfg.outputOrder)];
outDir = fileparts(cfg.outputFile);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end
writetable(T, cfg.outputFile, 'FileType', 'text', 'Delimiter', '\t');

fprintf('\nWrote %d unique samples: %s\n', height(T), cfg.outputFile);
fprintf('Feasible before deduplication: %d / %d (%.1f%%)\n', ...
        feasible, generated, 100*feasible/generated);
fprintf('Selection: %s; seed: %s; precision: %d decimals\n', ...
        cfg.selection, string(cfg.seed), cfg.decimalPlaces);
for j = 1:size(samples,2)
    v = samples(:,j);
    if max(v) == min(v)
        occupied = 1;
    else
        occupied = nnz(histcounts(v, linspace(min(v),max(v),11)));
    end
    fprintf('%s: observed [%.3g, %.3g], occupied deciles %d/10\n', ...
            cfg.outputOrder{j}, min(v), max(v), occupied);
end

%% LOCAL FUNCTIONS
function U = lhsUnit(n, d)
    U = zeros(n,d);
    for k = 1:d
        column = ((0:n-1)' + rand(n,1))/n;
        U(:,k) = column(randperm(n));
    end
end

function x = expandBound(x, n)
    if isscalar(x)
        x = repmat(x,n,1);
    else
        x = x(:);
    end
    assert(numel(x) == n, 'A bound or constraint returned the wrong row count.');
end

function s = roundStruct(s, names, digits)
    for k = 1:numel(names)
        s.(names{k}) = round(s.(names{k}), digits);
    end
end

function valid = checkBounds(s, cfg, n)
    valid = true(n,1);
    for k = 1:numel(cfg.sampleOrder)
        x = s.(cfg.sampleOrder{k});
        lo = expandBound(cfg.lower{k}(s,n), n);
        hi = expandBound(cfg.upper{k}(s,n), n);
        tol = 1e-9 * max(1, max(abs(lo),abs(hi)));
        valid = valid & isfinite(x) & isfinite(lo) & isfinite(hi) & ...
                x >= lo-tol & x <= hi+tol;
    end
end

function chosen = farthestPointSubset(pool, count, ranges)
    scaled = (pool-ranges(:,1)')./(ranges(:,2)-ranges(:,1))';
    m = size(pool,1);
    chosen = zeros(count,1);
    used = false(m,1);
    minDist2 = inf(m,1);
    chosen(1) = randi(m);
    for k = 1:count
        idx = chosen(k);
        used(idx) = true;
        minDist2 = min(minDist2, sum((scaled-scaled(idx,:)).^2,2));
        minDist2(used) = -inf;
        if k < count
            best = find(minDist2 == max(minDist2));
            chosen(k+1) = best(randi(numel(best)));
        end
    end
end

