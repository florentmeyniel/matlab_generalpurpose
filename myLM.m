function res = myLM(data, desmat)
% compute a multiple linear regression in parallel along several dimensions.
%
% Usage: res = myLM(data, desmat)
% data: matrix of data. Rows: replication, columns: data to be processed in parallel. 
% desmat: design matrix, rows: replication, column: regressors. 
%
% Florent Meyniel 2012-11-12

% remove nans
nanind_data = any(isnan(data), 2);
nanind_desm = any(isnan(desmat), 2);
nnind = ~(nanind_data | nanind_desm);

nn_data   = data(nnind, :);
nn_desmat = desmat(nnind, :);

% using code from regress.m
% Use the rank-revealing QR to remove dependent columns of X.
ndata        = size(nn_data, 2);
[n, nreg]    = size(nn_desmat);
[Q, R, perm] = qr(nn_desmat, 0);
p = sum(abs(diag(R)) > max(n, nreg)*eps(R(1)));
if p < nreg
    warning('stats:regress:RankDefDesignMat', ...
        'X is rank deficient to within machine precision.');
    R = R(1:p, 1:p);
    Q = Q(:, 1:p);
    perm = perm(1:p);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
b = zeros(nreg, ndata);
b(perm, :) = R \ (Q'*nn_data);

res = b; 