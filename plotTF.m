function plotTF(X, XTicks, XFlag, YTicks, YFlag, fig, title, ColorLim)
% plot the time-frequency plot, with lower frequency on the lower part of
% the graph
%
% Usage: plotTF(X, XTicks, XFlag, YTicks, YFlag, fig, title, ColorLim)
%
% X: [freq x time] data matrix. Freq increases over lines and time 
%    increases over columns
% XTicks / XFlag: to specify the X axis (optional)
% YTicks / YFlag: to specify the Y axis (optional)
% fig: figure handle (optional)
% title: window title (optional)
% ColorLim: [Low High] bound of the color sclare, or 'maxabs' (optional)

if ~exist('fig', 'var')
    fig = figure;
else
    if isempty(fig)
        fig = figure;
    else
        figure(fig)
    end
end

X = squeeze(X);
if ndims(X)>2
    error('even without singleton dimension, the data have %d > 2 dimensions', ndims(X))
end

set(gcf, 'Color', [1 1 1])

if exist('ColorLim', 'var') 
    if ischar(ColorLim) && strcmpi(ColorLim, 'maxabs')
        ColorLim = max(abs(X(:)))*[-1 1];
    end
    if size(X, 1) == numel(YFlag) && size(X, 2) == numel(XFlag)
        imagesc(YFlag, XFlag, X, sort(ColorLim))
    else
        imagesc(X, sort(ColorLim))
    end
else
    if size(X, 1) == numel(YFlag) && size(X, 2) == numel(XFlag)
        imagesc( X, sort(ColorLim))
    else
        imagesc(X)
    end
end
set(gca, 'YDir', 'Normal')

if exist('XTicks', 'var') && exist('XFlag', 'var')
    if length(XTicks) ~= length(XFlag)
        warning('XTicks and XFlag lengths differ')
    else
        
        % conver XFlag to a cell of strings
        % Strangely, with numeric values, when there is a different number 
        % of digits to print across ticks, numbers with fewer digits are
        % misaligned with the ticks.
        tmp = cell(1, numel(XFlag));
        for k = 1:numel(XFlag)
            tmp{k} = num2str(XFlag(k));
        end
        set(gca, 'XTick', XTicks,  'XTickLabel', tmp)
        % set(gca, 'XTick', XTicks,  'XTickLabel', XFlag)
        
    end
end

if exist('YTicks', 'var') && exist('YFlag', 'var')
    if length(YTicks) ~= length(YFlag)
        warning('YTicks and YFlag lengths differ')
    else
        set(gca, 'YTick', YTicks, 'YTickLabel', YFlag)
    end
end

if exist('title', 'var') && ischar(title)
   set(gcf, 'Name', title) 
end

xlabel('time')
ylabel('frequency')


