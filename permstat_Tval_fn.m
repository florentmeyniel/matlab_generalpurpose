function [t_thd, pval] = permstat_Tval_fn(x, alpha, np, tail, methodperm)
% Compute the t distribution under the null hypothesis using a permutation
% scheme. Report the t-value for which the Type 1 error for a more extreme 
% statistics is smaller than the alpha rate.
% 
% Usage:
%   [t_thd, pval] = permstat_Tval_fn(x, alpha, np, tail)
%
%     Input:
%           x: data matrix. The first dimension corresponds to replicates 
%              (e.g. subjects) and the other dimensions to the data for 
%              each replication
%       alpha: Type 1 error rate (default is 0.05)
%          np: amount of permutations (default: 10000). Set 'exact' to test
%          all permutations.
%        tail: 'right' (default), 'left' or 'both'. 'right' produce the 
%              test to give the probability of a higher t-value by chance. 
%              'both' give this probability for a bilateral test.  
%  methodperm: (optional) 'exact' to use exactely distinct permutations (but
%              time & memory demanding) or 'approx' to compute random
%              permutations (Default = approx)
%           
%    Output:
%       t_thd: minimum t-value controlled at the alpha rate
%       p_val: p-value corresponding to the t-value computed on the data
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
end
if ischar(np) && ~strcmp(np, 'exact')
    np = 10000;
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
% compute the closet t-value for a higher statistic at the alpha level
t_thd = prctile(nullt(:), alpha);

% compute p-values
if nargout > 1    
    N = np*ns;
    maxnullt = max(nullt(:));
    snullt = sort(nullt(:));
    
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

