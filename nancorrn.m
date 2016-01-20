function [c, p] = nancorrn(X)
% call corr(X with X a matrix, but after removing nans from X

nanind = all(~isnan(X), 2);
[c, p] = corr(X(nanind,:));
