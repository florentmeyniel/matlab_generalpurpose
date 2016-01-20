function [EffortDur, RestDur] = SplitRestEffort(data, ...
                                                method, ...
                                                ON2OFFDeadtime, ...
                                                OFF2ONDeadtime, ...
                                                thd, ...
                                                SR, ...
                                                sd_threshold)

                                            
% Function that implement the method of rest and effort duration.
% It works on the data of a trial, and the duration of each rest and effort
% period selected.

% INITIALIZATION
% ==============

EffortDur = [];
RestDur   = [];

% FIND ONSET AND OFFSET
% =====================

% Use the FindOnOff.m function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[On Off] = FindOnOff(data.grip, ...
    'ON2OFFDeadtime', ON2OFFDeadtime, ...
    'OFF2ONDeadtime', OFF2ONDeadtime, ...
    'thd', thd, ...
    'SR', SR, ...
    'sd_threshold', sd_threshold);

% Add ON and OFF so that each trial starts with an ON and ends with an OFF
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% case 1: trial starts with an ON and ends with an OFF
% => nothing to do
if ~isempty(On) ...                                                         % = there is an ON
        && ((~isempty(Off) && isempty(Off(Off<On(1)))) || isempty(Off)) ... % = ON is the first event
        && ~isempty(Off) ...                                                % = there is an OFF
        && ((~isempty(On) && isempty(On(On>Off(end)))) || isempty(On))      % = OFF is the last event
end

% case 2: trial starts with an ON and ends with an ON
% => put an OFF at the end of the trial
if ~isempty(On) ...                                                         % = there is an ON
        && ((~isempty(Off) && isempty(Off(Off<On(1)))) || isempty(Off)) ... % = ON is the first event
        && ~isempty(On) ...                                                 % = there is an ON in trial n
        && (isempty(Off) || (~isempty(Off) && isempty(Off(Off>On(end)))))   % = ON is the last event
    
    Off = [Off, length(data.grip)];
    
end

% case 3: trial starts with an OFF and ends with an OFF
% => put an ON at the start of the trial
if ~isempty(Off) ...                                                        % = there is an OFF
        && (isempty(On) || (~isempty(On) && isempty(On(On<Off(1)))))...     % OFF is the first event
        && ~isempty(Off) ...                                                % = there is an OFF
        && ((~isempty(On) && isempty(On(On>Off(end)))) || isempty(On))      % = OFF is the last event
    
    On = [1, On];
end

% case 4: trial starts with an OFF and ends with an ON
% => put an OFF at the end of the trial and a ON at the
% begining
if ~isempty(Off) ...                                                        % = there is an OFF
        && (isempty(On) || (~isempty(On) && isempty(On(On<Off(1)))))...     % OFF is the first event
        && ~isempty(On) ...                                                 % = there is an ON in trial n
        && (isempty(Off) || (~isempty(Off) && isempty(Off(Off>On(end)))))   % = ON is the last event
    
    On = [1, On];
    Off = [Off, length(data.grip)];
end

% convert on/offset to time (s), 0 is the trial
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if mean(data.time(end)-data.time(1)) > 1000           % time is in ms
    tOn  = (data.time(On)-data.time(1))/1000;
    tOff = (data.time(Off)-data.time(1))/1000;
else                                                  % time is in s
    tOn  = (data.time(On)-data.time(1));
    tOff = (data.time(Off)-data.time(1));
end

% SELECT EFFORT AND REST DEPENDING ON THE METHOD
% ==============================================
switch method
    
    case 'full'
        
        % -- Effort --
        for n_On=tOn
            EffortDur = [EffortDur, tOff(find(tOff(:)>n_On, 1, 'first')) - n_On];
        end
        
        % -- Rest --
        if tOn>0 % i.e. there is not an ON at trial onset
            RestDur = [RestDur, tOn(1)];
        end
        for n_Off=tOff
            if  ~isempty(tOn(tOn(:)>n_Off))
                RestDur = [RestDur, tOn(find(tOn(:)>n_Off, 1, 'first')) - n_Off];
            end
        end
        
    case 'min'
        
        % -- Effort --
        for n_On=tOn(1:end-1)
            EffortDur = [EffortDur, tOff(find(tOff(:)>n_On, 1, 'first')) - n_On];
        end
        if ~isempty(Off) && Off(end)~=length(data.grip) % i.e. last OFF was not added
            EffortDur = [EffortDur, tOff(end) - tOn(end)];
        end
        
        % -- Rest --
        for n_Off=tOff
            if  ~isempty(tOn(tOn(:)>n_Off))
                RestDur = [RestDur, tOn(find(tOn(:)>n_Off, 1, 'first')) - n_Off];
            end
        end
        
    case 'min+end'
        
        % -- Effort --
        for n_On=tOn(1:end)
            EffortDur = [EffortDur, tOff(find(tOff(:)>n_On, 1, 'first')) - n_On];
        end
        
        % -- Rest --
        for n_Off=tOff
            if  ~isempty(tOn(tOn(:)>n_Off))
                RestDur = [RestDur, tOn(find(tOn(:)>n_Off, 1, 'first')) - n_Off];
            end
        end
        
    case 'min+start'
        
        % -- Effort --
        for n_On=tOn(1:end-1)
            EffortDur = [EffortDur, tOff(find(tOff(:)>n_On, 1, 'first')) - n_On];
        end
        if ~isempty(Off) && Off(end)~=length(data.grip) % i.e. last OFF was not added
            EffortDur = [EffortDur, tOff(end) - tOn(end)];
        end
        
        % -- Rest --
        if tOn>0 % i.e. there is not an ON at trial onset
            RestDur = [RestDur, tOn(1)];
        end
        for n_Off=tOff
            if  ~isempty(tOn(tOn(:)>n_Off))
                RestDur = [RestDur, tOn(find(tOn(:)>n_Off, 1, 'first')) - n_Off];
            end
        end
        
    case 'min+start+coin'
        
        % -- Effort --
        for n_On=tOn(1:end-1)
            EffortDur=[EffortDur, tOff(find(tOff(:)>n_On, 1, 'first')) - n_On];
        end
        if ~isempty(Off) && Off(end)~=length(data.grip) % i.e. last OFF was not added
            EffortDur = [EffortDur, tOff(end) - tOn(end)];
        end
        
        % -- Rest --
        if tOn>0 % i.e. there is not an ON at trial onset
            % add one second to 1st rest
            RestDur = [RestDur, tOn(1) + 1]; %1: duration (s) of coin display
        end
        for n_Off=tOff
            if  ~isempty(tOn(tOn(:)>n_Off))
                RestDur = [RestDur, tOn(find(tOn(:)>n_Off, 1, 'first')) - n_Off];
            end
        end
        
    case 'full+coin'
        
        % -- Effort --
        for n_On=tOn
            EffortDur = [EffortDur, tOff(find(tOff(:)>n_On, 1, 'first')) - n_On];
        end
        
        % -- Rest --
        if tOn>0 % i.e. there is not an ON at trial onset
            RestDur = [RestDur, tOn(1) + 1]; %1: duration (s) of coin display
        end
        for n_Off=tOff
            if  ~isempty(tOn(tOn(:)>n_Off))
                RestDur = [RestDur, tOn(find(tOn(:)>n_Off, 1, 'first')) - n_Off];
            end
        end
end