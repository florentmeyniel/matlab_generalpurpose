function [t, p, df] = Tval_uqM_uqV(X1, X2)
% In short: same as ttest2 for unequal variance, but also return t-value &
% degree of freedom.
% 
% X1 & X2 are data (replicates x data) with unequal number of replicates 
% and unequal and variance.
% Compute the t-value(s) (corresponding to X1 - X2), the degrees of fredom 
% and the p-value for a more extreme statistics (bilateral test).
% 
% Usage
% [t, p, df] = Tval_uqM_uqV(X1, X2)
%          
%       t: t-value
%       p: p-value (bilateral test)
%      df: approximated degree of freedom
%
% NB: as df is approximated from the sample variances and sizes, df
% may be different for different data set even though they have the same
% sizes.
%
% Florent Meyniel 2012-10-23

% check size
if size(X1, 2) ~= size(X2, 2)
    error('X1 has %d and X2 has %d data point. SHOULD BE EQUAL!', size(X1, 2), size(X2, 2))
end

% compute mean
m1 = nanmean(X1, 1);
m2 = nanmean(X2, 1);

% compute variance
v1 = nanvar(X1, [], 1);
v2 = nanvar(X2, [], 1);

% compute replicate
n1 = size(X1, 1);
n2 = size(X2, 1);

% COMPUTE T-VALUE
t = (m1 - m2) ./ (sqrt( v1/n1 + v2/n2));

% COMPUTE DEGREE OF FREEDOM
% Welch-Satterthwaite equation to approximate the t-value significance
% (copmuted with unpooled variance) onto an ordinary Student's t
% distribution.
df = (v1/n1 + v2/n2).^2 ./ ( (v1/n1).^2/(n1-1) + (v2/n2).^2/(n2-1)  );

% COMPUTE SIGNIFICANCE
% this is the same as [~, p] = ttest2(X1, X2, [], [], 'unequal')
p = 2*tcdf(-abs(t), df);
 



