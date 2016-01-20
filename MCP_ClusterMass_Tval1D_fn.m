function [Cluster, StatAtFDR] = MCP_ClusterMass_Tval1D_fn(x, t_thd, FDR, nrand, method, methodperm)
% Non parametric estimation of cluster False Discovery Rate, testing for
% the data significance away from 0. This aims at addressing the multiple
% comparison issue, i.e of comparing the significance away from 0 across
% mutliple time bins. 
%
% NB: the algorithm use a permutation strategy that avoid any repetition.
%
% Usage: [Cluster, StatAtFDR] = MCP_ClusterMass_Tval1D_fn(x, t_thd, FDR, nrand, method, methodperm)
%
% Output
%    Cluster{i}.p_unilpos : 0< FDR <1 for a higher cluster stat (unilateral)
%    Cluster{i}.p_unilneg : 0< FDR <1 for a lower cluster stat (unilateral)
%    Cluster{i}.p_bilat   : 0< FDR <1 for a lower & higher cluster stat (bilateral)
%    Cluster{i}.ind       : indices of the cluster in the time series
%    Cluster{i}.t         : t-values of the cluster in the time series
%    Cluster{i}.stat      : cluster stat
% 
%    StatAtFDR(1, 1) : cluster stat at FDR for a higher stat (unilateral)
%    StatAtFDR(1, 2) : cluster stat at FDR for a lower stat (unilateral)
%    StatAtFDR(2, :) : cluster stat interval (lower, higher) at FDR for a more extreme stat (bilateral)
%
% Input
%    x      : [nreplicate x ntimes] matrix of time series. replicate can be subjects for instance
%    t_thd  : t-value threshold to define cluster on the t-values computed from x
%    FDR    : (optional) False Discovery Rate (0< FDR < 1), i.e. the fraction of
%             clusters (defined using t_thd) of having, by chance, a higher cluster
%             statistics (Default = 0.05)
%    nrand  : (optional) amount of randomization to estimate the null-hypothesis
%             distribution (Default = 10000)
%    method : (optional) cluster statistic. Can be 'size' (cluster size) or 'weight'
%             (Default, cluster sum of values)
% methodperm: (optional) 'exact' to use exactely distinct permutations (but
%             time & memory demanding) or 'approx' to compute random
%             permutations (Default = approx)
%
% Florent Meyniel 2012-01-30
% Modif. 2013-07-22: add methodperm.

% CHECK INPUTS
% ============
if ~exist('t_thd', 'var') ||  isempty(t_thd)
    error('threshold t-value should be defined for clustering')
else
    if t_thd<0
        t_thd = -t_thd;
    end
end
if ~exist('FDR', 'var') || isempty(FDR)
    FDR = 0.05;
else
    if FDR > 1 || FDR < 0
        error('FDR: %d!!! should be: 0 < FDR < 1', FDR)
    end
end
if ~exist('nrand', 'var') || isempty(nrand)
    nrand = 10000;
end

if ~exist('method', 'var') || isempty(method)
    ClusterStatFun = @sum;
else
    if strcmpi('size', method)
       ClusterStatFun = @length;
    elseif strcmpi('weight', method)
        ClusterStatFun = @sum;
    else
        ClusterStatFun = @sum;
        warning('use cluster mass')
    end
end

if ~exist('method', 'var')
    methodperm = 'approx';
elseif ~(strcmp(methodperm, 'approx') || strcmp(methodperm, 'exact'))
        error('unknow tail type %s, must be ''approx'', ''exact''', methodperm)
end


% INITIALIZE VARIABLES
% ====================
nr      = size(x, 1); % amount of replicate
ns      = size(x, 2); % amount of data point (per replicate)
MaxStat = zeros(1, nrand);
MinStat = zeros(1, nrand);
AbsStat = zeros(1, nrand);

% compute permutations
switch methodperm
    case 'exact'
        try
            % Get all possible permutations
            allposs = ff2n(nr);
            allposs(allposs==0) = -1;
            
            % select a random subset of permutations
            nrandall = size(allposs, 1);
            indall = randperm(nrandall);
            randtable(:,:) = allposs(indall(1:nrand),:)';
            
        catch
            warning('exact method exceedes capacity -> use approx. instead')
            
            randper = double(rand(nr, nrand) > 0.5);
            randper(randper(:) == 0) = -1;
            randtable(:, :) = randper;
        end
    case 'approx'
        randper = double(rand(nr, nrand) > 0.5);
        randper(randper(:) == 0) = -1;
        randtable(:, :) = randper;
end

% COMPUTE CLUSTER STATISTIC DISTRIBUTION OVER PERMUTATIONS
% =======================================================
for i_rand = 1:nrand
    
    % randomize data
    x_rand = x .* (randtable(:, i_rand)*ones(1, ns));
    
    % compute the t-values
    t_rand = mean(x_rand, 1) ./ stderror(x_rand, 1);
    
    % Distribution for a higher statitics
    [clusterstart, clusterend, ncluster] = clusterize(t_rand > t_thd);
    clustat = zeros(1, ncluster);
    for i_cluster = 1:ncluster
        clustat(i_cluster) = ClusterStatFun(t_rand(clusterstart(i_cluster) : clusterend(i_cluster)));
    end
    if ~isempty(clustat)
            MaxStat(i_rand) = max(clustat);
    end
    
    % Distribution for a lower statitics
    [clusterstart, clusterend, ncluster] = clusterize(t_rand < -t_thd);
    clustat = zeros(1, ncluster);
    for i_cluster = 1:ncluster
        clustat(i_cluster) = ClusterStatFun(t_rand(clusterstart(i_cluster) : clusterend(i_cluster)));
    end
    if ~isempty(clustat)
            MinStat(i_rand) = min(clustat);
    end
    
    % Distribution for a more extreme statitics
    [clusterstart, clusterend, ncluster] = clusterize(abs(t_rand) > t_thd);
    clustat = zeros(1, ncluster);
    for i_cluster = 1:ncluster
        clustat(i_cluster) = ClusterStatFun(t_rand(clusterstart(i_cluster) : clusterend(i_cluster)));
    end
    if ~isempty(clustat)
        if max(clustat) >= max(abs(clustat)) 
            AbsStat(i_rand) = max(clustat);
        else
            AbsStat(i_rand) = -max(abs(clustat));
        end
    end 
end  

% DEFINE CLUSTER STATISTICS AT THE FDR
% ====================================

StatAtFDR(1,1) = prctile(MaxStat, 100*(1-FDR));
StatAtFDR(1,2) = prctile(MinStat, 100*FDR);
StatAtFDR(2,:) = prctile(AbsStat, 100*[FDR/2 1-FDR/2]);

% COMPUTE THE ORIGINAL CLUSTER STATISTICS
% =======================================

% t-values for original data
t = mean(x, 1) ./ stderror(x, 1);

% get data cluster
[clusterstart, clusterend, ncluster] = clusterize(t > t_thd);
clustat = zeros(1, ncluster);
Cluster = cell(ncluster, 1);
for i_cluster = 1:ncluster
    % get cluster statistics
    clustat(i_cluster) = ClusterStatFun(t(clusterstart(i_cluster) : clusterend(i_cluster)));
   
    % compute FDR of the cluster statistics
    Cluster{i_cluster}.p_unilpos = ClusterStatFun(MaxStat>clustat(i_cluster))/nrand;
    Cluster{i_cluster}.p_unilneg = ClusterStatFun(MinStat<clustat(i_cluster))/nrand;
    Cluster{i_cluster}.p_bilat   = ClusterStatFun(abs(AbsStat)>clustat(i_cluster))/nrand;
    
    % report cluster characteristics
    Cluster{i_cluster}.ind       = clusterstart(i_cluster) : clusterend(i_cluster);
    Cluster{i_cluster}.t         = t(clusterstart(i_cluster) : clusterend(i_cluster));
    Cluster{i_cluster}.stat      = clustat(i_cluster);
end

%%%%%%%%%%%%%%%%%%%%%%%%SUBFUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clusterstart, clusterend, ncluster] = clusterize(data)
    
    % data = vector of 0 and 1 value. The function clusterize the
    % neighbouring equal value (i.e cluster of 1s and 0s)
    
    if size(data, 1) < size(data, 2)
       data = data'; 
    end
    data         = [0; data; 0];
    cluster      = diff(data);
    ncluster     = length(find(cluster))/2;
    clusterstart = find(cluster==1);
    clusterend   = find(cluster==-1)-1;
end

end