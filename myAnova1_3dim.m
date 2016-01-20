function [p, F, df1, df2] = myAnova1_3dim(x)
% Compute a one-way anova and factorize the computation over dependant
% variables
% 
% NB1: code based on anova1.m (stat toolbox)
% NB2: this is NOT a within subject anova (the model does not consider that
% a given subject s(j) produce the vector of obersavation o(j) over
% conditions
%
% Usage [p, F, df1, df2] = myAnova1_3dim(x)
% x: [replicate x condition x dep. var] matrix of data
% p: p-value for a higher F statitics
% F: F-statistic
% df1: condition degrees of freedom (= size(x,2) - 1)
% df2: error degrees of freedom (= #condition x (#replicate -1)
%                              i.e size(x,2)*(size(x,1)-1) )
%
% Florent Meyniel 2012-02-20

s   = size(x);
r   = s(1);                         % number of replicate
c   = s(2);                         % number of columns (=conditions)
n   = s(3);                         % number of dep. var.
xr  = x;
mu  = mean(mean(xr, 1), 2);
xr  = xr - repmat(mu, [r , c, 1]);  % center to improve accuracy
xm  = squeeze(mean(xr, 1));         % column means

df1 = c-1;                          % Column degrees of freedom
df2 = c*(r-1);                      % Error degrees of freedom

RSS = sum(r*(xm.*xm));              % Regression Sum of Squares
TSS = sum(sum(xr.*xr, 1), 2);       % Total Sum of Squares

% % alternative: include the mean (but should be zero...)
% gm  = mean(xm);
% rgm = repmat(gm, [c, 1]);
% RSS = sum(r*((xm - rgm).*(xm - rgm)));
% rrgm = repmat(reshape(gm, [1 1 n]), [r, c, 1]);
% TSS = sum(sum((xr - rrgm).*(xr - rrgm), 1), 2);  

TSS = squeeze(TSS)';                % Total Sum of Squares
SSE = TSS - RSS;                    % Error Sum of Squares

if (df2 > 0)
   mse = SSE/df2;                   % mean square error
else
   mse = NaN;
end

if (SSE~=0)
    F = (RSS/df1) ./ mse;           % F statistic
    p = 1-fcdf(F,df1,df2);          % Probability of F given equal means.
elseif (RSS==0)                     % Constant Matrix case.
    F = NaN;
    p = NaN;
else                                % Perfect fit case.
    F = Inf;
    p = 0;
end

