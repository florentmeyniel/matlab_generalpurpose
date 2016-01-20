function [c, p] = nancorr(X, Y)
% call corr(X, Y) with X & Y vectors of same length, but after removing
% nans from X & Y

X = X(:);
Y = Y(:);

nanind = ~isnan(X) & ~isnan(Y);

[c, p] = corr(X(nanind), Y(nanind));
