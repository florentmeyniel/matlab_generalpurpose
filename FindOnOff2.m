function [Times_On Times_Off]=FindOnOff2(data, varargin)
% 
% Find onset & offset that are not after the feedback 
% and exclude offset or onset that are too close from one another
%
% Usage: [n_ON n_OFF]=FindOnOff(data, ...)
%   data: vector of measurements
%   Times_On: timing effort onsets (i-th sample)
%   Times_Off: timing of effort offsets (i-th sample)
%   
%   Optional arguments:
%   'ON2OFFDeadtime': minimum dead time between an OFF and the previous ON
%       (s), default is 2s
%   'OFF2ONDeadtime': minimum dead time between an ON and the previous OFF
%       (s), dedault is 2s
%   'thd': threshold (ratio to Fmax, i.e. 1 is for Fmax) between rest and
%       effort, default is 0.5
%   'SR' : sampling rate (Hz), default is 1250Hz
%   'sd_threshold' : the criterion to find onset or offset, is that the derivate of
%       the signal is higher (absolute value) than the global standard
%       deviation * sd_threshold
%
%   CHANGE FROM VERSION 1
%   => if the trial ends with an effort period and that the maximal force
%   during this period is < tdh, the according Onset is discarded.
%
%
% Example:
%   [Times_On Times_Off]=FindOnOff(data, 'thd', 0.1)
% 
% Florent Meyniel

% Check input argument
if mod((nargin-1),2)~=0
    error('Wrong input arguments')
end
if nargin>1
    checkspell=0;
    for i = 1:length(varargin)
        checkspell = checkspell+strcmpi(varargin(i), 'ON2OFFDeadtime')...
            + strcmpi(varargin(i), 'OFF2ONDeadtime') + strcmpi(varargin(i), 'thd')... 
            + strcmpi(varargin(i), 'SR') + strcmpi(varargin(i), 'sd_threshold');
    end
    if (nargin-1)/2~=checkspell
        error('check spelling of input arguments')
    end
end

% Default Parameters
SR             = 1250;
ON2OFFDeadtime = 2*SR;
OFF2ONDeadtime = 2*SR;
thd            = 0.5; 
Times_On       = []; 
Times_Off      = [];

% Assign parameter to their input value (if any)
for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'SR'); SR=varargin{k+1};
        if ~isnumeric(varargin{k+1}); error('input must be numeric'); end        
    end
end
for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'ON2OFFDeadtime'); ON2OFFDeadtime=varargin{k+1}*SR;
        if ~isnumeric(varargin{k+1}); error('input must be numeric'); end        
    end
    if strcmpi(varargin{k}, 'OFF2ONDeadtime'); OFF2ONDeadtime=varargin{k+1}*SR;
        if ~isnumeric(varargin{k+1}); error('input must be numeric'); end        
    end
    if strcmpi(varargin{k}, 'thd'); thd=varargin{k+1};
        if ~isnumeric(varargin{k+1}); error('input must be numeric'); end        
    end
    if strcmpi(varargin{k}, 'sd_threshold'); sd_threshold=varargin{k+1};
        if ~isnumeric(varargin{k+1}); error('input must be numeric'); end        
    end
end

[Onset Offset]       = SearchFunction(data);
[Times_On Times_Off] = reject(Onset, Offset);
% assignin('base', 'Offset', Offset);
% assignin('base', 'Onset', Onset);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reject effort onset or offsets that are too close from one another
    function [Onset_kept, Offset_kept] = reject(Onset, Offset)
       % reject effort onset or offsets that are too close from one another 
        
        if isempty(Onset)
            Onset_kept  = []; 
            Offset_kept = [];
            return
        end
        
        Onset_kept  = Onset(1); 
        Offset_kept = [];
        
        %disp(['Onset ' num2str(length(Onset)) ' Offset ' num2str(length(Offset))])
        if (length(Offset)-length(Onset))>0
            disp(['Onset ' num2str(length(Onset)) ' Offset ' num2str(length(Offset))])
            error('odd number of offset')
        end
        
        if length(Onset)>length(Offset) % trial ends with a Onset marker
            for i = 2:length(Onset)
                if Onset(i)-Offset(i-1)>OFF2ONDeadtime
                    Onset_kept  = [Onset_kept Onset(i)];
                end
                if Offset(i-1)-Onset(i-1)>ON2OFFDeadtime
                    Offset_kept = [Offset_kept Offset(i-1)];
                end
            end
        else % trial ends with an Offset marker
            for i = 1:length(Offset)-1
                if Onset(i+1)-Offset(i)>OFF2ONDeadtime
                    Onset_kept  = [Onset_kept Onset(i+1)];
                end
            end
            for i = 1:length(Offset)
                if Offset(i)-Onset(i)>ON2OFFDeadtime
                    Offset_kept = [Offset_kept Offset(i)];
                end
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Onset, Offset] = SearchFunction(data)
        % find all effort onsets and offsets
        
        derivate           = diff(data);
        criterion          = std(diff(data))*sd_threshold;
        OnsetLib           = []; 
        OffsetLib          = [];
        Onset              = []; 
        Offset             = [];        
        lastonsetisliberal = 0;
        
        % 1- Define event liberaly
        % ~~~~~~~~~~~~~~~~~~~~~~~~
        % Look for effort onsets and offsets liberaly: some are actually
        % not offset or onset
        for i = 2:length(derivate)
            if abs(derivate(i))>criterion ...           % high speed
                    && abs(derivate(i-1))<criterion...  % transition from low speed
                    && derivate(i)<0 ...                % downward dynamics
                    && data(i)>thd                      % force level is high
                OffsetLib = [OffsetLib, i];             % => effort offset
            end
            if abs(derivate(i))>criterion ...           % high speed
                    && abs(derivate(i-1))<criterion...  % transition from low speed
                    && derivate(i)>0 ...                % upward dynamics
                    && data(i)<thd                      % force level is low
                OnsetLib = [OnsetLib, i];               % => effort onset
            end
        end
        
        if data(end)<thd % Add liberaly the last relax moment as an effort onset
            lastonsetisliberal = 1;
            OnsetLib = [OnsetLib, i];
        end
        
        % 2- Refine selection of Offsets
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Select the liberaly-defined offsets in a more conservative way: only
        % one offset is taken for each effort period (the last)
        for i = 1:length(OnsetLib)-1
            Offset = [Offset, ...
                max(OffsetLib(OffsetLib(:)>OnsetLib(i) & OffsetLib(:)<OnsetLib(i+1)))];
        end
        
        % 3- Refine selection of Onsets
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Select the liberaly-defined onsets in a more conservative way:
        % eliminate the last one if choosen liberaly and take only one onset
        % for each relax period
        if lastonsetisliberal==1
            OnsetLib = OnsetLib(1:end-1);
        end
        if ~isempty(OnsetLib)
            Onset = [Onset OnsetLib(1)];
        end
        
        % new code start: when many OnsetLib, choose the one for which the
        % force level is the lowest
        for i = 1:length(Offset)-1
            InterOff   = OnsetLib(OnsetLib(:)>Offset(i) & OnsetLib(:)<Offset(i+1));
            [~, ind]   = min(data(InterOff));
            Onset      = [Onset, InterOff(ind)];
        end
        if ~isempty(Offset) && ~isempty(OnsetLib(OnsetLib(:)>Offset(end)))
            AfterOff = OnsetLib(OnsetLib(:)>Offset(end));
            [~, ind] = min(data(AfterOff));
            Onset    = [Onset, AfterOff(ind)];
        end
        % new code end
        
        % old code: when many OnsetLib, choose the first one
%                 for i=1:length(OnsetLib)-1
%                     if ~isempty(Offset(Offset(:)>OnsetLib(i) & Offset(:)<OnsetLib(i+1)))
%                         Onset=[Onset OnsetLib(i+1)];
%                     end
%                 end
        % old code end
        
        % CHANGE SPECIFIC TO VERSION 2
        % if the last 'effort period' is on average < thd, it is not
        % further consider as an effort period.
        if ~isempty(Onset)
            if isempty(Offset) || isempty(find(Offset > Onset(end), 1))
                if max(data(Onset(end):end)) < thd
                    Onset = Onset(1:end-1);
                end
            end
        end
        
    end

end
