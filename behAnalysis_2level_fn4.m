function [bEi bRi bE bR] = behAnalysis_2level_fn4(EffortDur, RestDur, method, doplot, stat)

% function to compute a 2-level model analysis to estimate effects of 
% incentive, difficulty level, bloc, trial and within trial position on 
% effort and rest durations.
% At the 1st level, each subject is submitted to a multilinear regression.
% Beta coefficients are then put in a second-level, statistical analysis.
% Model with and without interaction are both estimated.
% The data are plotted and the regression coefficients are return
%
% Usage: [bEi bRi bE bR] = behAnalysis_2level_fn(EffortDur, RestDur, method, doplot, stat)
%   - EffortDur: matrix of effort duration (from an extractdata function)
%   - RestDur: matrix of rest duration (from an extractdata function)
%   - method: 'indiv' or 'group' method to plot data point of each subject or not
%   - doplot: 1 to do plot, 0 to avoid (default)
%   - bEi: regression coefficient for effort, with interactions
%   - bRi: regression coefficient for rest, with interactions
%   - bE: regression coefficient for effort, without interactions
%   - bR: regression coefficient for rest, without interactions
% 
% Regression coefficients: matrix regressor x subject
% Regressor order:
%   - Without interaction:
%       * constant
%       * RW
%       * FC
%       * Session
%       * Trial
%       * within trial
%   - With interactions:
%       * constant
%       * RW
%       * FC
%       * Session
%       * Trial
%       * within trial
%       * RW x Trial
%       * RW x Session
%       * RW x within trial
%       * FC x Trial
%       * FC x Session
%       * FC x within trial
%       * RW x FC
%
%
% variation wrt behAnalysis_2level_fn: 
% interactions are not orthogonolized, but estimated after main effect in
% their null space AND separately for each kind of fatigue.
% 
% Florent Meyniel 2011-09-01

% WITHOUT INTERACTION
% ===================

if ~exist('doplot', 'var'); doplot = 0; end
if ~exist('stat', 'var'); stat = 1; end
if doplot
    nfig = gcf;
end

for i_type=1:2
    clear p_val CI_val
    
    % define data
    if i_type==1
        data=[RestDur, ones(size(RestDur, 1),1)];
        if stat == 1
        fprintf('\n================= REST =================')
        end
    else
        data=[EffortDur, ones(size(EffortDur, 1),1)];
        if stat == 1
        fprintf('\n================ EFFORT ================')
        end
    end

    b_allSub=[];

    % --- FIRST LEVEL ---
    % ===================
    
    for i_sub=unique(data(:,1))'
        
        data_sub=data(data(:,1)==i_sub,:);
        
        % Add a column (8th) for period order within a trial
        % NB: within trial order is zscored trial-wise
        for i_block=1:max(data_sub(:,3))
            for i_hand=1:max(data_sub(:,2))
                for i_trial=1:max(data_sub(:,4))
                    data_sub(data_sub(:,2)==i_hand...
                             & data_sub(:,3)==i_block...
                             & data_sub(:,4)==i_trial...
                            ,8)=...
                        zscore(1:length(data_sub(data_sub(:,2)==i_hand...
                                                 & data_sub(:,3)==i_block...
                                                 & data_sub(:,4)==i_trial...
                                                 ,7)))';
                end
            end
        end

        % Zscore regressors by bloc & Hand
        X_bloc=zeros(size(data_sub,1),1);
        for i_hand=1:max(data_sub(:,2))
            X_bloc(data_sub(:,2)==i_hand)=...
                zscore(data(data_sub(:,2)==i_hand,3));
        end
        
        X_trial=zeros(size(data_sub,1),1);
        for i_hand=1:max(data_sub(:,2))
            for i_block=1:max(data_sub(:,3))
                X_trial(data_sub(:,3)==i_block & data_sub(:,2)==i_hand)=...
                    zscore(data_sub(data_sub(:,3)==i_block & data_sub(:,2)==i_hand,4));
            end
        end
        
        % Zscore regressors by bloc & Hand
        X_RW=zeros(size(data_sub,1),1);
        for i_hand=1:max(data_sub(:,2))
            for i_block=1:max(data_sub(:,3))
                X_RW(data_sub(:,3)==i_block & data_sub(:,2)==i_hand)=...
                    zscore(data_sub(data_sub(:,3)==i_block & data_sub(:,2)==i_hand,5));
            end
        end
        
        X_FC=zeros(size(data_sub,1),1);
        for i_hand=1:max(data_sub(:,2))
            for i_block=1:max(data_sub(:,3))
                X_FC(data_sub(:,3)==i_block & data_sub(:,2)==i_hand)=...
                    zscore(data_sub(data_sub(:,3)==i_block & data_sub(:,2)==i_hand,6));
            end
        end
        
        X_cst=ones(size(data_sub,1),1);
        X_intrial=data_sub(:,8);

        % multiple regression on native data
        [b]=regress(data_sub(:,7), ...
            [X_cst ...
            X_RW X_FC ...
            X_bloc X_trial X_intrial]);

        b_allSub=[b_allSub, b];
    end

    % --- SECOND LEVEL ---
    % ====================
    
    % comute ttest values on beta values
    for i_reg=2:size(b_allSub,1)
        [H P CI]=ttest(b_allSub(i_reg,:));
        p_val(i_reg-1)=P;
        CI_val(:,i_reg-1)=CI;
    end

    
    
    significance = cell(1,length(p_val));
    for i = 1:length(p_val)
        significance{i} = '   ';
        if p_val(i) < 0.001
            significance{i} = '***';
        elseif p_val(i) < 0.01
            significance{i} = '** ';
        elseif p_val(i) < 0.05
            significance{i} = ' * ';
        elseif p_val(i) < 0.1
            significance{i} = ' . ';
        end
    end
    if stat == 1
    % Sum up stats in the command line window
    fprintf('\n======= RESULTS NO INTERACTIONS ========\n\n')
    fprintf('          RW: p= %6.2e   %s\n', p_val(1), significance{1})
    fprintf('          FC: p= %6.2e   %s\n', p_val(2), significance{2})
    fprintf('       BLOCK: p= %6.2e   %s\n', p_val(3), significance{3})
    fprintf('       TRIAL: p= %6.2e   %s\n', p_val(4), significance{4})
    fprintf('WITHIN TRIAL: p= %6.2e   %s\n', p_val(5), significance{5})
    end
    % plot
    if doplot
        fig=figure(nfig + i_type);
        if i_type==1
            set(fig, 'Name', 'Mean rest duration')
        else
            set(fig, 'Name', 'Mean effort duration')
        end
        
        bar(nanmean(b_allSub(2:end,:),2), 'LineWidth', 2, 'FaceColor', 'b')
        hold on
        
        if strcmp(method, 'indiv')
            for i_sub=1:size(b_allSub(2:end,:),2)
                plot(1:length(b_allSub(2:end,i_sub)), b_allSub(2:end,i_sub), 'xg', 'MarkerSize', 10, 'LineWidth', 2)
            end
        end
        
        errorbar(1:size(nanmean(b_allSub(2:end,:),2),1), nanmean(b_allSub(2:end,:),2)', ...
            nanmean(b_allSub(2:end,:),2)'-CI_val(1,:), CI_val(2,:)'-nanmean(b_allSub(2:end,:),2), ...
            '.r', 'LineWidth', 3)
        set(gca, 'XTickLabel',{'RW','FC','BLOC','TRIAL','IN-TRIAL'},'XTick',[1 2 3 4 5], 'FontSize', 13, 'FontWeight', 'bold')
        set(gcf, 'Color', 'w')
    end
    
    if i_type==1
        bR = b_allSub;
    else
        bE = b_allSub;
    end

end


% WITH INTERACTIONS
% ===================
for i_type=1:2
    clear p_val CI_val
    
    % define data
    if i_type==1
        data=[RestDur, ones(size(RestDur, 1),1)];
        if stat == 1
        fprintf('\n================= REST =================')
        end
    else
        data=[EffortDur, ones(size(EffortDur, 1),1)];
        if stat == 1
        fprintf('\n================ EFFORT ================')
        end
    end

    b_allSub=[];

    % --- FIRST LEVEL ---
    % ===================
    
    for i_sub=unique(data(:,1))'
        
        data_sub=data(data(:,1)==i_sub,:);
        
        % Add a column (8th) for period order within a trial
        % NB: within trial order is zscored trial-wise
        for i_block=1:max(data_sub(:,3))
            for i_hand=1:max(data_sub(:,2))
                for i_trial=1:max(data_sub(:,4))
                    data_sub(data_sub(:,2)==i_hand...
                             & data_sub(:,3)==i_block...
                             & data_sub(:,4)==i_trial...
                            ,8)=...
                        zscore(1:length(data_sub(data_sub(:,2)==i_hand...
                                                 & data_sub(:,3)==i_block...
                                                 & data_sub(:,4)==i_trial...
                                                 ,7)))';
                end
            end
        end

        % Zscore regressors by bloc & Hand
        X_bloc=zeros(size(data_sub,1),1);
        for i_hand=1:max(data_sub(:,2))
            X_bloc(data_sub(:,2)==i_hand)=...
                zscore(data(data_sub(:,2)==i_hand,3));
        end
        
        X_trial=zeros(size(data_sub,1),1);
        for i_hand=1:max(data_sub(:,2))
            for i_block=1:max(data_sub(:,3))
                X_trial(data_sub(:,3)==i_block & data_sub(:,2)==i_hand)=...
                    zscore(data_sub(data_sub(:,3)==i_block & data_sub(:,2)==i_hand,4));
            end
        end
        
       % Zscore regressors by bloc & Hand
        X_RW=zeros(size(data_sub,1),1);
        for i_hand=1:max(data_sub(:,2))
            for i_block=1:max(data_sub(:,3))
                X_RW(data_sub(:,3)==i_block & data_sub(:,2)==i_hand)=...
                    zscore(data_sub(data_sub(:,3)==i_block & data_sub(:,2)==i_hand,5));
            end
        end
        
        X_FC=zeros(size(data_sub,1),1);
        for i_hand=1:max(data_sub(:,2))
            for i_block=1:max(data_sub(:,3))
                X_FC(data_sub(:,3)==i_block & data_sub(:,2)==i_hand)=...
                    zscore(data_sub(data_sub(:,3)==i_block & data_sub(:,2)==i_hand,6));
            end
        end
        
        X_cst=ones(size(data_sub, 1),1);
        X_intrial=data_sub(:,8);
        
        % include interactions
        X_RWxtrial=zscore((X_RW + abs(min(X_RW)) + 1).*X_trial);
        
        
        
        X_RWxbloc     = zscore((X_RW + abs(min(X_RW)) + 1).*X_bloc);
        X_RWxintrial  = zscore((X_RW + abs(min(X_RW)) + 1).*X_intrial);
        X_FCxtrial    = zscore((X_FC + abs(min(X_FC)) + 1).*X_trial);
        X_FCxbloc     = zscore((X_FC + abs(min(X_FC)) + 1).*X_bloc);
        X_FCxintrial  = zscore((X_FC + abs(min(X_FC)) + 1).*X_intrial);
        X_RWxFC       = zscore((X_RW + abs(min(X_RW)) + 1).*X_FC);

        % multiple regression on native data
        b = zeros(13, 1);
        desmat1 = [X_cst ...
            X_RW X_FC...
            X_bloc X_trial X_intrial ...
            ];
        
        desmatRW = [X_RWxtrial X_RWxbloc X_RWxintrial];
        
        desmatFC = [X_FCxtrial X_FCxbloc X_FCxintrial];
        
        desmatRWFC = [X_RWxFC];
        
        [b(1:6), ~, R] = regress(data_sub(:,7), desmat1);
        
        [b(7:9)] = regress(R, desmatRW);
        [b(10:12)] = regress(R, desmatFC);
        [b(13)] = regress(R, desmatRWFC);

        b_allSub=[b_allSub, b];
    end

    % --- SECOND LEVEL ---
    % ====================
    
    % comute ttest values on beta values
    for i_reg=2:size(b_allSub,1)
        [H P CI]=ttest(b_allSub(i_reg,:));
        p_val(i_reg-1)=P;
        CI_val(:,i_reg-1)=CI;
    end

    

    significance = cell(1,length(p_val));
    for i = 1:length(p_val)
        significance{i} = '   ';
        if p_val(i) < 0.001
            significance{i} = '***';
        elseif p_val(i) < 0.01
            significance{i} = '** ';
        elseif p_val(i) < 0.05
            significance{i} = ' * ';
        elseif p_val(i) < 0.1
            significance{i} = ' . ';
        end
    end
    
    if stat == 1
    % Sum up stats in the command line window
    fprintf('\n======= RESULTS WITH INTERATIONS =======\n\n')
    fprintf('             RW: p= %6.2e   %s\n', p_val(1), significance{1})
    fprintf('             FC: p= %6.2e   %s\n', p_val(2), significance{2})
    fprintf('          BLOCK: p= %6.2e   %s\n', p_val(3), significance{3})
    fprintf('          TRIAL: p= %6.2e   %s\n', p_val(4), significance{4})
    fprintf('   WITHIN TRIAL: p= %6.2e   %s\n', p_val(5), significance{5})
    fprintf('       RW*Trial: p= %6.2e   %s\n', p_val(6), significance{6})
    fprintf('        RW*Bloc: p= %6.2e   %s\n', p_val(7), significance{7})
    fprintf('RW*Within Trial: p= %6.2e   %s\n', p_val(8), significance{8})
    fprintf('       FC*Trial: p= %6.2e   %s\n', p_val(9), significance{9})
    fprintf('        FC*Bloc: p= %6.2e   %s\n', p_val(10), significance{10})
    fprintf('FC*Within Trial: p= %6.2e   %s\n', p_val(11), significance{11})
    fprintf('          RW*FC: p= %6.2e   %s\n', p_val(12), significance{12})
    end
    
    % plot
    if doplot
        fig=figure(nfig + 2 + i_type);
        if i_type==1
            set(fig, 'Name', 'Mean rest duration, with interactions')
        else
            set(fig, 'Name', 'Mean effort duration, with interactions')
        end
        bar(nanmean(b_allSub(2:end,:),2), 'LineWidth', 2, 'FaceColor', 'b')
        hold on
        
        if strcmp(method, 'indiv')
            for i_sub=1:size(b_allSub(2:end,:),2)
                plot(1:length(b_allSub(2:end,i_sub)), b_allSub(2:end,i_sub), 'xg', 'MarkerSize', 10, 'LineWidth', 2)
            end
        end
        
        % plot(b_allSub(2:end,:), '.g')
        errorbar(1:size(nanmean(b_allSub(2:end,:),2),1), nanmean(b_allSub(2:end,:),2)', ...
            nanmean(b_allSub(2:end,:),2)'-CI_val(1,:), CI_val(2,:)'-nanmean(b_allSub(2:end,:),2), ...
            '.r', 'LineWidth', 3)
        set(gca, 'XTickLabel',{'RW','FC','Bl','TR','ITr', 'RW*Tr', 'RW*Bl', 'RW*iTr',...
            'FC*Tr', 'FC*Bl', 'FC*iTr', 'RW*FC'},'XTick',[1:12], 'FontSize', 13, 'FontWeight', 'bold')
        set(gcf, 'Color', 'w')
    end
    
    if i_type==1
        bRi = b_allSub;
    else
        bEi = b_allSub;
    end

end

