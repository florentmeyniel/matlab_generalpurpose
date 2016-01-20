function plotfigbeha_fn(EffortDur, RestDur, axisE, axisR, nYTicks)
% Compute and plot the intersubject mean and standard error of effort and rest 
% durations for each factor level
% 
% Usage: plotfigbeha_fn(EffortDur, RestDur, axisE, axisR, nYTicks)
%   - EffortDur: matrix of effort duration (from an extractdata function)
%   - RestDur: matrix of rest duration (from an extractdata function)
%   - axisE: plot axis limit for Effort ([] sets it to Matlab defaults)
%   - axisR: plot axis limit for Rest ([] sets it to Matlab defaults)
%   - nYTicks: number of ticks labels ([] sets it to 3)
% 
% Florent Meyniel 2011-09-01

if isempty(nYTicks)
    nYTicks = 3;
end

nsub = length(unique(EffortDur(:,1)));

% COMPUTE INTERSUBJECT MEAN AND VARIANCE
% ======================================

% 1- Subject Specific Mean
MrFc = zeros(nsub, 3); MrRw = zeros(nsub, 3);
MeFc = zeros(nsub, 3); MeRw = zeros(nsub, 3);
for irw=1:3
    for ifc=1:3
        for isub=1:nsub
            
            MeFc(isub, ifc) = nanmean(...
                EffortDur(EffortDur(:,1)==isub & ...
                              EffortDur(:,6)==ifc...
                              , 7));
            
            MeRw(isub, irw) = nanmean(...
                EffortDur(EffortDur(:,1)==isub & ...
                              EffortDur(:,5)==irw...
                              , 7));
                
            MrFc(isub, ifc) = nanmean(...
                RestDur(RestDur(:,1)==isub & ...
                            RestDur(:,6)==ifc...
                            , 7));
            
            MrRw(isub, irw) = nanmean(...
                RestDur(RestDur(:,1)==isub & ...
                            RestDur(:,5)==irw...
                            , 7));            
        end
    end
end

% 2- mean and SE over subjects

mMrFc = mean(MrFc, 1);
mMeFc = mean(MeFc, 1);
mMrRw = mean(MrRw, 1);
mMeRw = mean(MeRw, 1);

eMrFc = stderror(MrFc, 1);
eMeFc = stderror(MeFc, 1);
eMrRw = stderror(MrRw, 1);
eMeRw = stderror(MeRw, 1);

% PLOT DATA
% =========

fig=figure;
set(fig, 'Name', 'Inter Subject')
set(fig, 'Color', 'w')

subplot(2,2,1) % Effort Cost
bar(mMeFc, 'LineWidth', 1.5, 'FaceColor', 'w', 'EdgeColor', 'r'), hold on, 
errorbar(mMeFc, eMeFc, 'r.', 'LineWidth', 1), ylabel('Effort Duration (s)', 'FontSize', 12, 'FontWeight', 'normal')
set(gca, 'XTick', [1:3], 'XTickLabel', {'70%', '80%', '90%'}, 'FontSize', 12, 'FontWeight', 'normal')
if ~isempty(axisE)
    axis(axisE)
    set(gca, 'YTick', axisE(4).*linspace(0, 1, nYTicks), 'YTickLabel', cellstr(num2str(axisE(4).*linspace(0, 1, nYTicks)'))', 'FontSize', 12, 'FontWeight', 'normal')
end

subplot(2,2,2) % Effort Benefit
bar(mMeRw, 'LineWidth', 1.5, 'FaceColor', 'w', 'EdgeColor', 'r'), hold on, 
errorbar(mMeRw, eMeRw, 'r.', 'LineWidth', 1)
set(gca, 'XTick', [1:3], 'XTickLabel', {'10c', '20c', '50c'}, 'FontSize', 12, 'FontWeight', 'normal')
if ~isempty(axisE)
    axis(axisE)
    set(gca, 'YTick', axisE(4).*linspace(0, 1, nYTicks), 'YTickLabel', cellstr(num2str(axisE(4).*linspace(0, 1, nYTicks)'))', 'FontSize', 12, 'FontWeight', 'normal')
end

subplot(2,2,3) % Rest Cost
bar(mMrFc, 'LineWidth', 1.5, 'FaceColor', 'w', 'EdgeColor', 'b'), hold on
errorbar(mMrFc, eMrFc, 'b.', 'LineWidth', 1),  ylabel('Rest Duration (s)', 'FontSize', 12, 'FontWeight', 'normal')
set(gca, 'XTick', [1:3], 'XTickLabel', {'70%', '80%', '90%'}, 'FontSize', 12, 'FontWeight', 'normal')
if ~isempty(axisR)
    axis(axisR)
    set(gca, 'YTick', axisR(4).*linspace(0, 1, nYTicks), 'YTickLabel', cellstr(num2str(axisR(4).*linspace(0, 1, nYTicks)'))', 'FontSize', 12, 'FontWeight', 'normal')
end

subplot(2,2,4) % Rest Benefit
bar(mMrRw, 'LineWidth', 1.5, 'FaceColor', 'w', 'EdgeColor', 'b'), hold on
errorbar(mMrRw, eMrRw, 'b.', 'LineWidth', 1)
set(gca, 'XTick', [1:3], 'XTickLabel', {'10c', '20c', '50c'}, 'FontSize', 12, 'FontWeight', 'normal')
if ~isempty(axisR)
    axis(axisR)
    set(gca, 'YTick', axisR(4).*linspace(0, 1, nYTicks), 'YTickLabel', cellstr(num2str(axisR(4).*linspace(0, 1, nYTicks)'))', 'FontSize', 12, 'FontWeight', 'normal')
end
end