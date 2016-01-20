function CheckSamplingOfRegressors(SPM, iOns, iP, TR)
% Plot an fMRI convolved regressors at the microtime resolution and at the
% MRI sampling rate.
% Note that the regressors are convolved by the HRF, but were not
% orthogonolize (they are as provided by the user, and store in 
% SPM.Sess.U.u
%
% Usage: 
%   CheckSamplingOfRegressors(SPM, iOns, iP, TR)
%       - SPM: an SPM structure (in which the design matrix was computed)
%       - iOns: the onset regressor number
%       - iP: the parametric modulation of this regressor (iP=1 is the
%         onset regressor itself. 
%       - TR (optional): a repetition time (for comparision purpose with
%         the TR used in the SPM structure).

% Sparse input
if nargin == 3
    TR = [];
end
    
% COMPUTE REGRESSOR & SAMPLING
% ============================

% Get regressors at microtime resolution
X = SPM.Sess.U(iOns).u;
X = full(X(:, iP));

% compute convolved regressor
BOLD = conv(X, SPM.xBF.bf)';

% Get the number of data point per scan from the SPM structure
BinPerScan1 = SPM.xBF.T;

% compute with another sampling if requiered
if ~isempty(TR)
    BinPerScan2 = TR/SPM.xBF.dt;
end

scanvec1 = round(1:BinPerScan1:length(X));

if ~isempty(TR)
    scanvec2 = round(1:BinPerScan2:length(X));
end

% Get the time vector
timevec = SPM.xBF.dt*(1:length(X));

% PLOT RESULT
% ===========

figure; set(gcf, 'Color', [1 1 1])

plot(timevec, BOLD(1:length(timevec)), 'LineWidth', 2)

hold on
plot(timevec(scanvec1), BOLD(scanvec1), 'o')
if ~isempty(TR)
    plot(timevec(scanvec2), BOLD(scanvec2), 'ro')
    legendtext = {'microtime res.', ...
        sprintf('chosen TR: %4.3f', SPM.xBF.T*SPM.xBF.dt), ...
        sprintf('other TR: %4.3f', TR)};
else
    legendtext = {'microtime res.', ...
        sprintf('chosen TR: %4.3f', SPM.xBF.T*SPM.xBF.dt)};
end

legend(legendtext, 'Location', 'Best')

plot(timevec(scanvec1), BOLD(scanvec1))
if ~isempty(TR)
    plot(timevec(scanvec2), BOLD(scanvec2), 'r')
end

ylabel('Predicted BOLD')
xlabel('time (s)')
if iP == 1
    % get name of onset regressor
    title(SPM.Sess.U(iOns).name(iP))
else
    % get name of parametric modulation
    title(SPM.Sess.U(iOns).P(iP-1).name)
end

