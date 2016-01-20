function [CI] = StudentCI(data, M, alpha)
%
% CI = StudentCI(data, M, alpha)
%  returns a 100*(1-ALPHA)% confidence interval for the true mean of X
%  being M. 
% defaut: M = 0
%         alpha = 0.05
%
% Florent Meyniel

if ~exist('M', 'var');     M     = 0    ; end
if ~exist('alpha', 'var'); alpha = 0.05 ; end

[~, ~, CI] = ttest(data, M, alpha);


end