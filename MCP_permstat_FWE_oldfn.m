function [t_FWE, pval] = MCP_permstat_FWE_fn(x, alpha, np)
% Compute the t distribution under the null hypothesis using a permutation
% scheme. Report the t-value for which the familywise error rate for a 
% higher statistic (i.e. one-side t-test) is controlled at the alpha level.
%
% The FWE is computed using the maximum statistics distribution (i.e. the
% distribution of the maximum t-value from each permutation).
% The maximul statistics control weakly and strongly for FWE. 
% See p 423 in 
%   Thomas Nichols & Satoru Hayasaka
%    Controlling the familywise error rate in functional neuroimaing: A
%    comparative review
%    Statistical Methods in Medical Research 2003; 12: 419-446
% 
% Usage:
%   [t_FWE, pval] = MCP_permstat_FWE_fn(x, alpha, np)
%     Input:
%           x: data matrix. The first dimension corresponds to replicates 
%              (e.g. subjects) and the other dimensions to the data for 
%              each replication
%       alpha: FWE rate (default is 0.05)
%          np: amount of permutations (default: 10000)
%           
%    Output:
%       t_FWE: minimum t-value controlled at FWE alpha rate
%       p_val: FWE p-value corresponding to the t-value computed on the data
%
% Florent Meyniel 2012-02-09


% CHECK INPUTS
% ============
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 95;
else
    if alpha > 1 || alpha < 0
        error('alpha: %d!!! should be: 0 < alpha < 1', alpha)
    else
        alpha = 100 * (1-alpha);
    end
end
if ~exist('np', 'var') || isempty(np)
    np = 10000;
end

% RESHAPE DATA
% ============
nr      = size(x, 1);  % amount of replicate
ns      = numel(x)/nr; % amount of data point (per replicate)
% collapse replicate over dimensions
if ndims(x) > 2
    x   = reshape(x, [nr, ns])'; 
end

% COMPUTE THE NULL DISTRIBUTION
% =============================
% Get all possible permutations
allposs = ff2n(nr);
allposs(allposs==0) = -1;

% select a random subset of permutations
nrandall = size(allposs, 1);
indall = randperm(nrandall);
randtable(:,1,:) = allposs(indall(1:np),:)';

% compute HO, the randomized data
r_r = repmat(randtable, [1 ns 1]);
r_x = repmat(x, [1, 1, np]);
H0   = r_r.*r_x;

% comute null t-distribution

nullt = squeeze(mean(H0, 1) ./ stderror(H0, 1));

% REPORT THE STATS
% ================
MaxTH0 = max(nullt);
t_FWE  = prctile(MaxTH0, alpha);

% compute p-values
if nargout > 1
    N      = np;
    maxnullt = max(MaxTH0);
    snullt = sort(MaxTH0(:));
    
    % t-value of the data
    t = mean(x, 1) ./ stderror(x, 1);
    
    count    = zeros(size(t));
    for it = 1:numel(t)
        if maxnullt > t(it) 
        count(it) = find(snullt > t(it), 1, 'first');
        else
            count(it) = N;
        end
    end    
    pval = 1 - count/N;    
end

