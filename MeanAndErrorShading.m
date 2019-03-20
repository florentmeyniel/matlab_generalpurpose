function [h] = MeanAndErrorShading(xval, M, SEM, alpha, colm, colsd, LWm, LWsd, contourstyle)
% Plot M and +/- SEM as a shaded area.
%
% Usage: MeanAndErrorShading(xval, M, SEM, alpha, colm, colsd, LWm, LWsd, contourstyle)
%   xval - the [1 x n] vector of mean values (it can be [])
%   M - the [1 x n] vector of mean values
%   SEM - the [1 x n] vector of error values
%   alpha - transparenct coefficient of the shadding (>0 and <1)
%   colm - the color [R G B] for the mean
%   colsd - the color [R G B] of the shading ([] set to colm)
%   LWm - the line width of the mean ([] set to 1)
%   LWsd - the line width of the mean ([] set to 1)
%   contourstyle - the contour style '-', '--', ':', '-.' or 'none' ([] set
%   to 'none')
% 
%   h is the "plot" handle (useful for figure legend when this function is
%   called multiple times.
%
% Florent Meyniel 2011-08-16

if isempty(colsd)
    colsd = colm;
end
if isempty(LWm)
    LWm = 1;
end
if isempty(LWsd)
    LWsd = 1;
end
if isempty(contourstyle)
    contourstyle = 'none';
end

if size(M,1)>size(M,2)
    M = M';
end
if size(SEM,1)>size(SEM,2)
    SEM = SEM';
end

AZ = M + SEM;
BZ = M - SEM;

AZ_ind = find(~isnan(AZ));
AZ = AZ(AZ_ind);
BZ = BZ(AZ_ind);

if ~isempty(xval)
    AZ_val = xval(AZ_ind);
else
    AZ_val = AZ_ind; 
end

fill([AZ_val flipud(AZ_val')'],[AZ flipud(BZ')'],...
    'k', ...                            % what is this line for? but it doesn't work without it...
    'EdgeColor', colsd, ... 
    'LineWidth',LWsd,...
    'LineStyle',contourstyle,...
    'FaceColor',colsd,...
    'FaceAlpha',alpha);
hold on
h = plot(AZ_val, M, 'LineWidth', LWm, 'Color', colm);

