function Pval = permstat_Tval2Pval_fn(x, Tval, np, methodperm)
% Compute the t distribution under the null hypothesis using a permutation
% scheme. Report the p-value corresponding to the input t-value. 
% 
% Usage:
%   [t_thd, pval] = permstat_Tval_fn(x, Tval, np)
%
%     Input:
%           x: data matrix. The first dimension corresponds to replicates 
%              (e.g. subjects) and the other dimensions to the data for 
%              each replication
%        Tval: T-value
%          np: (optional) amount of permutations (default: 10000)
%  methodperm: (optional) 'exact' to use exactely distinct permutations (but
%              time & memory demanding) or 'approx' to compute random
%              permutations (Default = approx)
%           
%    Output:
%       Pval: p-value corresponding to the T-value
%       
% Florent Meyniel 2012-02-09
% Modif. 2013-07-22: add methodperm.


% CHECK INPUTS
% ============
if ~exist('np', 'var') || isempty(np)
    np = 10000;
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

% REPORT THE STATS
% ================
snullt = sort(nullt(:));
N = np*ns;
if isempty(find(snullt > Tval, 1, 'first'))
    Pval = 0;
else
    Pval = 1 - find(snullt > Tval, 1, 'first')/N;
end

