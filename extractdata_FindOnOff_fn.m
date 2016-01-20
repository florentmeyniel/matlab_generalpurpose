function [EffortDur, RestDur, FB, PERF]=extractdata_FindOnOff_fn(opt)
%
% This function loads data from subjects, block, hands and trials specified,
% and compute all effort and rest durations.
% Effort onset (ON) and offset  (OFF) are computed using the FindOnOff
% function.
%
% USAGE
% -----
% [EffortDur, RestDur, FB, PERF] = extractdata_FindOnOff_fn(opt)
%
% OPTIONS TO DEFINE
% ----------------
% 
%        opt.datadir : full path of data directory (string)
%     opt.FilePrefix : prefix of .mat file in the data directory (string)
%      opt.subIDlist : subject number list (vector) or cell of strings
%     opt.def_blocks : set block number (vector) or 'default' (to take all)
%      opt.def_hands : set hands (e.g. ['d', 'g'])
%     opt.def_trials : set trials (vector) or 'default' (to take all)
% opt.ON2OFFDeadtime : minimum dead time between an OFF and the previous ON (s)
% opt.OFF2ONDeadtime : minimum dead time between an ON and the previous OFF (s)
%            opt.thd : threshold (ratio to Fmax, i.e. 1 is for Fmax) between rest and effort
%             opt.SR : sampling rate (Hz)
%   opt.sd_threshold : the criterion to find onset or offset, is that the derivate of
%                       the signal is higher (absolute value) than the global standard
%                       deviation * sd_threshold
%         opt.method: 'full', 'min', 'min+end', 'min+start' , 'min+start+coin' or 'full+coin'
%                       5 methods are proposed:
%                                 'min' : the initial rest period and the last effort period when
%                                         truncated by trial end are rejected
%                             'min+end' : same as 'min' but the last effort period is not rejected
%                                         when truncated
%                           'min+start' : same as 'min', but the intial rest period is not
%                                         rejected
%                                'full' : same as 'min+start' and 'min+end' altogether
%                           'full+coin' : same as full but include coin display in the 1st rest 
%                       'min+start+coin : same as 'min+start' but include coin display in the 1st rest
%
% OUTPUT DATA 
% -----------
%
% columns:
%   1- subject number
%   2- hand (1: opt.def_hand(1), 2: opt.def_hand(2))
%   3- block number
%   4- trial number
%   5- reward level (1: lowest, 3:highest)
%   6- force level (1: lowest, 3:highest)
% (6'- Diplayed force level, if different from applied) 
%   7- data, can be:
%       * duration (ms) for EffortDur and RestDur
%       * FB (€) 
%       * time step for PERF
%
% NB: suppose that data file name is: 
% FilePrefix, subIDlist, def_hands, 'bloc', block#, '.mat'
%
% 
% function version of extractdata_FindOnOff
% Florent Meyniel 2011-07-04


datadir         = opt.datadir;
FilePrefix      = opt.FilePrefix;
subIDlist       = opt.subIDlist; 
def_blocks      = opt.def_blocks;
def_hands       = opt.def_hands;
def_trials      = opt.def_trials;
ON2OFFDeadtime  = opt.ON2OFFDeadtime;
OFF2ONDeadtime  = opt.OFF2ONDeadtime;
thd             = opt.thd;
SR              = opt.SR;
sd_threshold    = opt.sd_threshold;
method          = opt.method;

% check
methodlist = {'min', 'min+end', 'min+start', 'full', 'full+coin', 'min+start+coin'};
if ~strcmp(method, methodlist)
    error('The specified medhod ''%s'' is not recognized. It should be one of: \n ''%s'' \n ''%s'' \n ''%s'' \n ''%s'' \n ''%s'' \n ''%s''' ,...
        method, methodlist{1}, methodlist{2}, methodlist{3}, methodlist{4}, methodlist{5}, methodlist{6})
end

if ~iscell(subIDlist)
    newsubIDlist = cell(1, length(subIDlist));
    for iSub = 1:length(subIDlist)
        newsubIDlist{iSub} = num2str(subIDlist(iSub));
    end
    subIDlist = newsubIDlist;
end
    
% --- COMPUTE DURATIONS ---
% =========================

EffortDur = [];
RestDur   = []; 
FB        = []; 
PERF      = [];

for sub = 1:length(subIDlist)
    for hand = 1:length(def_hands)
        
        if strcmp(def_blocks, 'default')
            n = size(dir(strcat(datadir, FilePrefix, subIDlist{sub}, def_hands(hand), '*.mat')),1);
            blocks_list = 1:n;
           if n == 0
               fprintf('\nWARNING: could not find files for subject #%s \n', subIDlist{sub});
           end
        else
            blocks_list = def_blocks;
        end
            
        
        for block = blocks_list
            
            data = load(strcat(datadir, FilePrefix, ...
                subIDlist{sub}, def_hands(hand), 'bloc', num2str(block), '.mat'));
            
            if isfield(data, 'gripdata1')
                data.gripdata = data.gripdata1;
            elseif isfield(data, 'gripdata2')
                data.gripdata = data.gripdata2;
            end
            
            if isfield(data, 'data1')
                data.data = data.data1;
            elseif isfield(data, 'data2')
                data.data = data.data2;
            end
                
            if strcmp(def_trials, 'default')
                trial_list = 1: size(data.gripdata, 2);
            else
                trial_list = def_trial;
            end
            
            for ntrial = trial_list
                
                [trialEffortDur, trialRestDur] = ...
                    SplitRestEffort(data.gripdata{ntrial}, ...
                    method, ...
                    ON2OFFDeadtime, ...
                    OFF2ONDeadtime, ...
                    thd, ...
                    SR, ...
                    sd_threshold);
                
                trialCond = [sub, hand, block, ntrial, data.data(ntrial, 4:end-2)];
                
                if ~isempty(trialEffortDur)
                    EffortDur = [EffortDur; repmat(trialCond, length(trialEffortDur), 1), trialEffortDur'];
                end
                if ~isempty(trialRestDur)
                    RestDur   = [RestDur;   repmat(trialCond, length(trialRestDur), 1),   trialRestDur'];
                end
                
                PERF = [PERF; trialCond, data.data(ntrial, end-1)];
                FB   = [FB;   trialCond, data.data(ntrial, end)];
                
            end
        end
    end
end

end