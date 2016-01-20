function [EffortDur, RestDur] = extractdata_markers_fn(opt)

%
% This function loads .ds file data for the specified subjects.
% Effort onset (ON) and offset  (OFF) are not computed in this function,
% but simply set to the EMG_ON and EMG_OFF markers.
%
% USAGE
% -----
% [EffortDur, RestDur] = extractdata_marker_fn(opt)
%
% OPTIONS TO DEFINE
% ----------------
% - opt.H.Fs: sampling frequency
% - opt.H.nSamplesPre: number of pre-trigger time point (i.e. before t=0)
% - opt.datadir: full path of data directory (string) containing subject
%       folders
% - opt.FilePrefix: prefix of .ds files in the data directory (string)
% - opt.subIDlist: subject number list (vector) 
% - opt.progressbar: =1 to display the progress bar, =0 not to display it
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
%   7- data, can be duration (s) for EffortDur and RestDur
% 
% Florent Meyniel 2011-08-31

Cond = {'FC1RW1' 'FC2RW1' 'FC3RW1' 'FC1RW2' 'FC2RW2' 'FC3RW2' 'FC1RW3' 'FC2RW3' 'FC3RW3'};
code = [      1        1        1        2        2        2        3        3        3 ;... % RW
              1        2        3        1        2        3        1        2        3]';   % FC

% GET FULL PATH OF FILES
% ======================

subdir = dir(strcat(opt.datadir, 'sujet*'));
EffortDur = [];
RestDur   = [];

% chekc input
if length(opt.subIDlist) > length(subdir)
    error('%d subject directory and %d subject requested', length(opt.subIDlist), length(subdir))
end

if opt.progressbar
    pg_nStep = length(opt.subIDlist)*8;
    pg_Progress = 0;
    pg_h = waitbar(pg_Progress/pg_nStep, 'Computing durations...');
    set(pg_h, 'Color', [1 1 1])
end

for iSub = opt.subIDlist
    
    sessdir = dir(strcat(opt.datadir, '/', subdir(iSub).name, '/', opt.FilePrefix, '*'));
    
    for iSess = 1:length(sessdir)
        
        % Update progress bar
        if opt.progressbar
            pg_Progress = pg_Progress + 1;
            waitbar(pg_Progress/pg_nStep, pg_h)
        end
        
        % GET MARKER TIMINGS
        % ~~~~~~~~~~~~~~~~~~
        
        % set filename
        filename = strcat(opt.datadir, '/', subdir(iSub).name, '/', sessdir(iSess).name);
        
        % load timings
        Event = read_event_ds(filename);
        [TOn, ~]   = getMarkerTiming('EMG_ON', Event, opt.H);
        [TOff, ~]  = getMarkerTiming('EMG_OFF', Event, opt.H);
        [T_BL, ~]  = getMarkerTiming('BASELINE', Event, opt.H);
        [T_FB, ~]  = getMarkerTiming('FEEDBACK', Event, opt.H);
        
        CondTiming = zeros(length(Cond), 3);
        for iCond = 1:length(Cond)
            [CondTiming(iCond,3), ~] = getMarkerTiming(Cond(iCond), Event, opt.H);
            CondTiming(iCond, 1:2) = code(iCond, :);
        end
        [~, ind] = sort(CondTiming(:,3));
        CondTiming = CondTiming(ind, :);
        
        % get which hand was used
        if isempty(find(filename(end-3)==['2' '4' '6' '8'], 1));
            hand = 1; % i.e. right hand
        else hand = 0;
        end
        
        % get the session number
        nBloc = round(iSess/2);
        
        % COMPUTE EFFORT PERIODS
        % ~~~~~~~~~~~~~~~~~~~~~~
        
        % preallocate matrices
        dataEffort = zeros(length(TOn), 7);
        
        for iOn = 1:length(TOn)
            nextT_FB = T_FB(find(T_FB > TOn(iOn), 1, 'first'));
            prevT_BL = T_BL(find(T_BL < TOn(iOn), 1, 'last'));
            
            try
                prevT_BL(1); % check that there is a previous BL (cf. subject 7 run01 & 02)
                nextTOff = TOff(TOff > TOn(iOn) & TOff < (nextT_FB + 1)); % include a margin of 1s after FB
            catch
                nextTOff = [];
            end
            
            if ~isempty(nextTOff)
                nextTOff = nextTOff(1); % takes the first TOff after TOn
                
                % get trial number
                nTrial = length(find(T_BL < TOn(iOn)));
                                
                % get condition
                ind = find(CondTiming(:,3) < TOn(iOn), 1, 'last');
                
                % Report all results
                dataEffort(iOn, 1) = iSub;                    % # subject
                dataEffort(iOn, 2) = hand;                    % hand(1: right, 2: left)
                dataEffort(iOn, 3) = nBloc;                   % # bloc
                dataEffort(iOn, 4) = nTrial;                  % # Trial
                dataEffort(iOn, 5) = CondTiming(ind, 1);      % RW
                dataEffort(iOn, 6) = CondTiming(ind, 2);      % FC 
                dataEffort(iOn, 7) = nextTOff - TOn(iOn);     % duration
            end
        end    
        
        % COMPUTE REST PERIODS
        % ~~~~~~~~~~~~~~~~~~~~
        
        % preallocate matrice
        dataRest = zeros(length(TOn), 7);
        
        for iOn = 1:length(TOn)
            prevTOff = TOff(find(TOff < TOn(iOn), 1, 'last'));
            prevCondTiming = CondTiming(find(CondTiming(:,3) < TOn(iOn), 1, 'last'), 3);
            prevT_BL = T_BL(find(T_BL < TOn(iOn), 1, 'last'));
            
            if ~isempty(prevT_BL) % for subject 7: no ON/OFF excepted on BL onset in run 01 & 02. This if makes sure that nothing is done for these runs
                if(isempty(prevTOff) || prevTOff < prevCondTiming) % i.e. iOn is the trial 1st event
                    % set epoch bounds
                    duration = TOn(iOn) - prevCondTiming;
                    
                else  % i.e. iOn follows an Off belonging to the same trial
                    % set epoch bounds
                    duration = TOn(iOn) - prevTOff;
                end
                
                % get trial number
                nTrial = length(find(T_BL < TOn(iOn)));
                
                % get condition
                ind = find(CondTiming(:,3) < TOn(iOn), 1, 'last');
                
                % reports condition
                dataRest(iOn, 1) = iSub;                    % # subject
                dataRest(iOn, 2) = hand;                    % hand(1: right, 2: left)
                dataRest(iOn, 3) = nBloc;                   % # bloc
                dataRest(iOn, 4) = nTrial;                  % # Trial
                dataRest(iOn, 5) = CondTiming(ind, 1);      % RW
                dataRest(iOn, 6) = CondTiming(ind, 2);      % FC 
                dataRest(iOn, 7) = duration;                % duration
            end
        end
        
        % Concatenate data
        EffortDur = [EffortDur; dataEffort];
        RestDur   = [RestDur; dataRest];
        
    end
end

if opt.progressbar
    close(pg_h)
end

ind = EffortDur(:,1)~=0;
EffortDur = EffortDur(ind, :);
ind = RestDur(:,1)~=0;
RestDur = RestDur(ind, :);

end