function [t_FWE, pval] = MCP_permstat_FWE_fn(x, alpha, np, tail, methodperm)
% Compute the t distribution under the null hypothesis using a permutation
% scheme. Report the t-value for which the Type 1 error.
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
%   [t_FWE, pval] = MCP_permstat_FWE_fn(x, alpha, np, tail)
%     Input:
%           x: data matrix. The first dimension corresponds to replicates 
%              (e.g. subjects) and the other dimensions to the data for 
%              each replication
%       alpha: (optional) FWE rate (default is 0.05)
%          np: (optional) amount of permutations (default: 10000). Set 'all' to test
%              all permutations
%        tail: (optional) 'right' (default), 'left' or 'both'. 'right' produce the 
%              test to give the probability of a higher t-value by chance. 
%              'both' give this probability for a bilateral test.  
%  methodperm: (optional) 'exact' to use exactely distinct permutations (but
%              time & memory demanding) or 'approx' to compute random
%              permutations (Default = approx)
%    Output:
%       t_FWE: minimum t-value controlled at FWE alpha rate
%       p_val: FWE p-value corresponding to the t-value computed on the data
%
% Florent Meyniel 2012-02-09
% Modif. 2013-07-22: add methodperm.


% CHECK INPUTS
% ============
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 95;
else
    for k = 1:length(alpha)
        if alpha(k) > 1 || alpha(k) < 0
            error('alpha: %d!!! should be: 0 < alpha < 1', alpha(k))
        else
            alpha(k) = 100 * (1-alpha(k));
        end
    end
end
if ~exist('np', 'var') || isempty(np)
    np = 10000;
else
    if strcmpi(np, 'all')
        np = 2^size(x, 1);
    elseif isnumeric(np)
    else
        error('number of permutation not recognized')
    end
end

if ~exist('tail', 'var')
    tail = 'right';
else
    if ~(strcmp(tail, 'right') || strcmp(tail, 'left') || strcmp(tail, 'both'))
        error('unknow tail type %s, must be ''both'', ''right'' or ''left''', tail)
    end
end

if ~exist('methodperm', 'var')
    methodperm = 'approx';
elseif ~(strcmp(methodperm, 'approx') || strcmp(methodperm, 'exact'))
        error('unknow tail type %s, must be ''approx'', ''exact''', methodperm)
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
switch methodperm
    case 'exact'
        try
            % Get all possible permutations
            allposs = ff2n(nr);
            allposs(allposs==0) = -1;
            
            % select a random subset of permutations
            nrandall = size(allposs, 1);
            indall = randperm(nrandall);
            randtable(:,1,:) = allposs(indall(1:np),:)';
            
        catch
            warning('exact method exceedes capacity -> use approx. instead')
            
            randper = double(rand(nr, np) > 0.5);
            randper(randper(:) == 0) = -1;
            randtable(:, 1, :) = randper;
        end
    case 'approx'
        randper = double(rand(nr, np) > 0.5);
        randper(randper(:) == 0) = -1;
        randtable(:, 1, :) = randper;
end
    
% compute HO, the randomized data
r_r = repmat(randtable, [1 ns 1]);
r_x = repmat(x, [1, 1, np]);
H0   = r_r.*r_x;

% comute null t-distribution

nullt = squeeze(mean(H0, 1) ./ stderror(H0, 1));
if strcmp(tail, 'both')
    nullt = abs(nullt); % rectify distribution to achieve bilaterality
end
if strcmp(tail, 'left')
    nullt = -nullt; % reverse distribution to reach left tail    
end

% REPORT THE STATS
% ================
MaxTH0 = nanmax(nullt);
t_FWE  = prctile(MaxTH0, alpha);
if strcmp(tail, 'left')
    t_FWE = -t_FWE; % reverse distribution to reach left tail
end

% compute p-values
if nargout > 1
    N      = np;
    maxnullt = nanmax(MaxTH0);
    snullt = sort(MaxTH0(:));
    
    % t-value of the data
    t = nanmean(x, 1) ./ stderror(x, 1);
    if strcmp(tail, 'left')
        t = -t; % reverse distribution to reach left tail
    end
    if strcmp(tail, 'both')
        t = abs(t); % reverse distribution to reach left tail
    end
    
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

