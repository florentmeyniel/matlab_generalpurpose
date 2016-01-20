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
% ColorLim: [Low High] bound of the color sclae (optional)

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

if exist('ColorLim', 'var') && length(ColorLim) == 2
    imagesc(flipud(X), sort(ColorLim))
else
    imagesc(flipud(X))
end

if exist('XTicks', 'var') && exist('XFlag', 'var')
    if length(XTicks) ~= length(XFlag)
        warning('XTicks and XFlag lengths differ')
    else
        set(gca, 'XTick', XTicks,  'XTickLabel', XFlag)
    end
end


if exist('YTicks', 'var') && exist('YFlag', 'var')
    if length(YTicks) ~= length(YFlag)
        warning('YTicks and YFlag lengths differ')
    else
        set(gca, 'YTick', YTicks, 'YTickLabel', fliplr(YFlag))
    end
end

if exist('title', 'var') && ischar(title)
   set(gcf, 'Name', title) 
end

xlabel('time')
ylabel('frequency')


