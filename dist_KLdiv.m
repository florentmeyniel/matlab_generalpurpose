function KL = dist_KLdiv(Q, P, doWarning)
% Compute the Kullback–Leibler divergence between probability distributions 
% Q and P
% 
% Usage:
% KL = dist_KLdiv(Q, P, doWarning)
% NB: doWarning: 1 allow a warning message on whether the distribution is 
% normalized to sum to 1 or the divergence not defined (default: 0) 

if nargin == 2
    doWarning = 0;
end

% Check that both probability distribution sum to 1
if abs((sum(Q) - 1)) > 2*eps
    msg = sprintf(['the 1st distribution is not a probability distribution! '...
            '\n ... it is now normalized.']);
    Q = Q / sum(Q);
end
if abs((sum(P) - 1)) > 2*eps
    msg = sprintf(['the 2nd distribution is not a probability distribution! '...
        '\n ... it is now normalized.']);
    P = P / sum(P);
end

% check absolute continuity (when P is 0, then Q is 0)
indQ0 = find(Q == 0);
indP0 = find(P == 0);
if length(intersect(indP0, indQ0)) ~= length(indQ0)
    msg = sprintf(['KL not defined because at least one value satifies:', ...
        'P(i) = 0 & Q(i) ~= 0']);
    KL = NaN;
else
    ind = Q > 0;
    KL = sum(P(ind) .* (log(P(ind)) - log(Q(ind))), 2);
end

if doWarning == 1
    dist(msg)
end
