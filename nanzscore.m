function Z = nanzscore(X, varargin)
% 
% compute zscore of X (even with NaN values)
% 
% Usage: Z = stderror(X, dim)
%   - dim is the dimension used for computation
%      by default: dim=1 (over rows)
%      
% Florent Meyniel

% default option
dim=1;

% if X has only one dimension
if size(X,1)==1 || size(X,2)==1
    if size(X,1)>size(X,2)
        dim=1;
    else
        dim=2;
    end
end

% check varargin
if ~isempty(varargin)
    if size(varargin)>1
        error('only one optional argument is requiered')
    elseif ~isnumeric(varargin{1})
        error('dimension DIM must be numeric')
    else
        dim=varargin{1};
    end
    if dim>length(size(X))
        error('DIM is %d whereas the data dimensionality is %d', dim, length(size(X)))
    end
    if rem(dim,1)~=0
        error('DIM must be an natural number')
    end
end
if dim == 1 
    Z= (X - ones(size(X,1), 1)*nanmean(X, dim))./(ones(size(X,1), 1)*nanstd(X,[],dim));
else
    Z= (X - nanmean(X, dim)*ones(1, size(X,2)))./(nanstd(X,[],dim)*ones(1, size(X,2)));
end
end