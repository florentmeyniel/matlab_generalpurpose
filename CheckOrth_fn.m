function CheckOrth_fn(SPM, nSess, BF, preproc)
% Function equivalent of CheckOrth
%
% CheckOrth_fn(SPM, nSess, BF, preproc)
%   * SPM: the SPM structure (computed with spm_fMRI_design.m)
%   * nSess: the session number
%   * BF: the tempora basis function number.
%   * preproc: optional argument. 'w' (default) for the full preprocessing
%     of the design matrix, 'o' to stop just after orthogonization.
% Compare visually regressors as specified
% by the user and as used by SPM after orthogonalisation, whitening &
% filtering.
%

% get current figure number
FigNb = gcf;
try FigNb = FigNb.Number; end

if nargin < 4
    preproc = 'w';
else
    if ~ismember(preproc, {'w', 'o'})
        error('Unrecognized preprocessing mode %s. It should be ''o'' or ''w''', preproc)
    end
end

% check nSess
NbOfSess = size(SPM.Sess, 2);
if nSess>NbOfSess
    error('session rank can''t be > %d (= # of sessions)', NbOfSess)
end

Lscan   = SPM.nscan;            % number of scans per session
dt      = SPM.xBF.dt;           % microtime resolution in sec
dt2scan = SPM.xBF.T;            % # of microtime resolution data point in a scan
TR      = SPM.xBF.dt*SPM.xBF.T; % TR in sec
Lsec    = Lscan*TR;             % plot data until this time (sec)
nBF     = SPM.xBF.order;        % number of temporal basis functions

% get design matrix
if strcmp(preproc, 'w')
    desmat  = SPM.xX.xKXs.X;    % whitened design matrix
elseif strcmp(preproc, 'o')
    desmat  = SPM.xX.X;         % before whitening
end

% check BF
if BF>nBF
    error('basis function number can''t be > %d (= # of basis function)', nBF)
end

RegPerSess = []; % number of regressor per session
for iSess = 1:nSess
    RegPerSess(iSess) = length(SPM.Sess(iSess).col);
end

% column index of the regressors in the design matrix
RegRank = sum(RegPerSess(1:nSess-1))...
        + BF:nBF:sum(RegPerSess(1:nSess));

% scan range; i.e. row indices in the design matrix
ScanRange = sum(Lscan(1:nSess-1))+1 : sum(Lscan(1:nSess));

% estimate number of "dummy" scans
X = desmat(ScanRange, sum(RegPerSess(1:nSess-1))+ 1 : sum(RegPerSess(1:nSess)));
mini = [];
for iReg = 1:RegPerSess(nSess)-6; % 6 movement parameters
    mini(iReg) = find(X(:,iReg), 1, 'first');
end
dummyscans = min(mini);

nreg = 0;
for iU=1:size(SPM.Sess(1,nSess).U,2)
    
    % U REGRESSORS
    % ============
    
    nreg = nreg+1;
    
    % construct Reg at dt resolution using onset and durations
    reg_t = dt:dt:Lsec(nSess);
    reg_val = zeros(1,length(reg_t));
    for i = 1:length(SPM.Sess(1, nSess).U(1, iU).ons)
        if SPM.Sess(1, nSess).U(1, iU).dur(i)==0
            reg_val(reg_t>=SPM.Sess(1, nSess).U(1, iU).ons(i)...
                & reg_t<SPM.Sess(1, nSess).U(1, iU).ons(i)+dt)...
                = 1;
        else
            reg_val(reg_t>SPM.Sess(1, nSess).U(1, iU).ons(i) & ...
                reg_t<(SPM.Sess(1, nSess).U(1, iU).ons(i)+SPM.Sess(1, nSess).U(1, iU).dur(i)))...
                = 1;
        end
    end
    
    % Convolve with temporal basis function and resample at scan resolution
    reg_conv = conv(SPM.xBF.bf(:, BF), reg_val);
    reg_conv = reg_conv(1:length(reg_val));
    reg_conv_ds = reg_conv(1:dt2scan:end);
    
    % amplitude scaling
    scale = mean(abs(desmat(ScanRange, RegRank(nreg))))...
          /mean(abs(reg_conv_ds));
    
    % plot
    figure(FigNb+nreg)
    set(gcf, 'Name', SPM.Sess(1, nSess).U(1, iU).name{1})
    set(gcf, 'Color', 'w')
    plot(reg_conv_ds*scale, 'LineWidth', 2)
    hold on
    plot(desmat(ScanRange, RegRank(nreg)), 'r', 'LineWidth', 2)
    legend('designed', 'used'), xlabel('scan #')
    
    % PARAMETRIC MODULATIONS
    % ======================
    if size(SPM.Sess(1,nSess).U(1,iU).P,2)>1 || ~strcmp('none', SPM.Sess(1,nSess).U(1,iU).P.name)
        for iP = 1:size(SPM.Sess(1, nSess).U(1, iU).P,2)
            
            nreg = nreg+1;
            
            % construct Reg at dt
            reg_t = dt:dt:Lsec(nSess);
            reg_val = zeros(1, length(reg_t));
            for i = 1:length(SPM.Sess(1, nSess).U(1, iU).ons)
                pmod = zscore(SPM.Sess(1, nSess).U(1, iU).P(iP).P);
                if SPM.Sess(1,nSess).U(1,iU).dur(i) == 0
                    reg_val(reg_t>=SPM.Sess(1, nSess).U(1, iU).ons(i)...
                        & reg_t<SPM.Sess(1, nSess).U(1, iU).ons(i)+dt)...
                        = pmod(i);
                else
                    reg_val(reg_t>SPM.Sess(1, nSess).U(1, iU).ons(i) & ...
                        reg_t<(SPM.Sess(1, nSess).U(1, iU).ons(i)+SPM.Sess(1, nSess).U(1, iU).dur(i)))...
                        = pmod(i);
                end
            end
            
            % Convolve with temporal basis function and resample at scan resolution
            reg_conv = conv(SPM.xBF.bf(:, BF), reg_val);
            reg_conv((dummyscans-1)*dt2scan+1:end) = zscore(reg_conv((dummyscans-1)*dt2scan+1:end)); % get first dummy scans data out of scaling
            reg_conv = reg_conv(1:length(reg_val));
            reg_conv_ds = reg_conv(1:dt2scan:end);
            
            % amplitude scaling
            scale = mean(abs(desmat(ScanRange, RegRank(nreg))))...
                  /mean(abs(reg_conv_ds));
            
            % plot
            figure(FigNb+nreg)
            set(gcf, 'Name', strcat([SPM.Sess(1, nSess).U(1, iU).P(iP).name, ' (pmod)']))
            set(gcf, 'Color', 'w')
            plot(reg_conv_ds*scale, 'LineWidth', 2)
            hold on,
            plot(desmat(ScanRange, RegRank(nreg)), 'r', 'LineWidth', 2)
            legend('designed', 'used'), xlabel('scan #')
        end
    end
end
