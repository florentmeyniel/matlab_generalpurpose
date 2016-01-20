function [RR, PR, QR, RS, QS, RT...
          pkR, pkP, pkQ, pkS, pkT...
          Ramp] = GetCardiacComp(ECG, Fs, verbose)
% Compute the cardiac component timing: RR, PR, QR, RS, QS and RT segments
% and compute the R peak amplitude.
% The timing is in second; the amplitude in the unit of the recording.
%
% Note: the ECG is band-passed filter between 1Hz and 60Hz before
% examining the peaks. 
%
% Method:
%   1- compute a 1st-pass estimate of the R-peak, using a threshold-based
%   peak detection
%   2- use the 1st-pass estimate to epoch the data around the R-peaks and
%   estimate the canonical PQRST complex using PCA
%   3- use the (flipped and zero-mean) canonical PQRST complex as a kernel
%   to the convolve with the data
%   4- the peak of the convolution product is shifted to correspond to the
%   R-peak and the peaks are detected to using again, a threshold-based 
%   detection.
%   5- Q and S peaks are defined as minima around the R peak
%   6- P peaks are defined using template match for the P segment
%
% Usage: 
%     [RR, PR, QR, RS, QS, RT, ...
%           pkR, pkP, pkQ, pkS, pkT...
%          Ramp]= GetCardiacComp(ECG, Fs, verbose)
%
%   Inputs:
%         ECG: time series of the ECG
%          Fs: sampling rate of the time series
%     verbose: get the graph of PQRS complex from PCA, and ask whether to 
%              switch to mean when the PCA / mean correlation is low.
% 
%   Outputs:
% RR, PR, QR, RS, QS: duration of the segment (s)
% pkR, pkP, pkQ, pkS: sample index
%               Ramp: amplitude (recording unit)
%
% NB: RR = diff(pkR)/Fs; PR = (pkR - pkP)/Fs ...
%
%   Florent Meyniel 2012-02-07

if ~exist('verbose', 'var')
    verbose = 0;
end

% 1st PASS: ECG SIGNAL AND PEAK-BASED DECTETION
% =============================================

% FILTER DATA
hpf   = 1;                            % high pass frequency
[b a] = butter(4, hpf/(Fs/2), 'high');
ECG_f = filtfilt(b, a , ECG)'; 
lpf   = 60;                           % low pass frequency
[b a] = butter(4, lpf/(Fs/2), 'low');
ECG_f = filtfilt(b, a , ECG_f)'; 

% CHECK ORIENTATION
if median(ECG_f)>mean(ECG_f)
    ECG_f = -ECG_f;                     % check that R peak is on top
    warning('CHECK ORIENTATION!')
end

% FIND R PEAKS
pk = FindECGPeaks(ECG_f, Fs);           % peak position (in sample)

% 2nd PASS: CARDIAC CYCLE DETERMINE USING PCA
% ===========================================

% COMPUTE PCA TO FIND THE CANONICAL CARDIAC COMPLEX
% Epoch the data around the 1st-pass based R-peaks
befR   = 0.25;
aftR   = 0.4;
Lbef   = round(befR*Fs);
Laft   = round(aftR*Fs);
L      = Lbef+Laft+1;
Nstart = find((pk-Lbef)>0, 1, 'first');
Nend   = find(pk+Laft<=length(ECG_f), 1, 'last');
EPO    = zeros(L, length(pk(Nstart:Nend)));
for iEpoch = Nstart:Nend;
   EPO(:, iEpoch) = ECG_f(pk(iEpoch)-Lbef:pk(iEpoch)+Laft); 
end

% compute PCA
[P, T] = princomp(EPO);

% check the orientation
if mean(P(:,1))<0
    T(:,1) = -T(:,1);
end
template = T(:,1);

if corr(template, mean(EPO, 2)) < 0.9
    warning('low correlation 1st PCA component / mean = %4.3f', corr(template, mean(EPO, 2)))
end

if verbose == 1
    if corr(template, mean(EPO, 2)) >= 0.9
    plot([1:length(template)]/Fs, template)
    else
        subplot(1, 3, 1)
        plot([1:length(template)]/Fs, template) , title('PCA');
        subplot(1, 3, 2)
        plot([1:length(template)]/Fs, -template) , title('neg PCA');
        subplot(1, 3, 3)
        plot([1:length(template)]/Fs, mean(EPO, 2)) , title('mean');
        
        drawnow
        % Ask to choose
        button = questdlg('look figure!', 'Make your choice', 'PCA', 'negPCA', 'mean', 'PCA');
        if strcmp(button, 'negPCA')
            template = -template;
        end
        if strcmp(button, 'mean')
            template = mean(EPO, 2);
        end
        if strcmp(button, 'negmean')
            template = mean(EPO, 2);
        end        
    end
        
end
    
% SEARCH THE P, Q, R, S PEAKS
% ===========================

% R PEAKS
% -------
basisR  = flipud(template);               % reverse kernel to time direction
convolR = conv(ECG_f, basisR);            % convolve data and kernel
[~, basis_indmaxR] = max(basisR);         % get kernel maxima (= the R peak)
Rpeak_conv = convolR(basis_indmaxR:end);  % temporally realign to R-peak timing
pkR = FindECGPeaks(Rpeak_conv, Fs)';      % R peak position (in sample)

% Q AND S PEAKS
% --------------
QRlag = 0.05;
SRlag = QRlag;
pkQ   = zeros(1,length(pkR));
pkS   = zeros(1,length(pkR));

% check that there is enough signal for the 1st peak
if pkR(1)- round(QRlag *Fs) < 0
    pkP(1) = NaN;
    firstiR = 2;
else firstiR = 1;
end

% search peak using the extrema around the R peak on the ECG
for iR = firstiR:length(pkR)
    % Q peak
    [~, tmp] = min(ECG_f(pkR(iR)-round(QRlag*Fs) : pkR(iR)));    
    pkQ(iR)  = tmp + pkR(iR)-round(QRlag *Fs);
    
    % R peak
    [~, tmp] = min(ECG_f(pkR(iR) : pkR(iR)+round(SRlag*Fs)));    
    pkS(iR) = tmp + pkR(iR);
end

% P PEAKS
% -------
% search peak using template match and relative position to R peak
PRlag  = 0.05;
basisP = T(1:round((befR - PRlag)*Fs),1);
basisP = flipud(basisP);                  % reverse kernel to time direction
convolP = conv(ECG_f, basisP);            % convolve data and kernel
[~, basis_indmaxP] = max(basisP);         % get kernel maxima (= the P peak)
Ppeak_conv = convolP(basis_indmaxP:end);  % temporally realign to P-peak timing

% check that there is enough signal for the 1st peak
pkP = zeros(1,length(pkR));
if pkR(1) - round(befR *Fs) < 0
    pkP(1) = NaN;
    firstiR = 2;
else firstiR = 1;
end

% search peak using the extrema around the R peak on the ECG
for iR = firstiR:length(pkR)
    [~, tmp] = max(Ppeak_conv(pkR(iR)-round(befR *Fs) : pkR(iR)-round(PRlag*Fs)));
    pkP(iR) = tmp + pkR(iR)-round(befR *Fs);
end

% T PEAKS
% -------
% search peak using template match and relative position to R peak
RTlag  = 0.1;
basisT = T(round((befR + RTlag)*Fs:end),1);
basisT = flipud(basisT);                  % reverse kernel to time direction
convolT = conv(ECG_f, basisT);            % convolve data and kernel
[~, basis_indmaxT] = max(basisT);         % get kernel maxima (= the T peak)
Tpeak_conv = convolT(basis_indmaxT:end);  % temporally realign to T-peak timing

% search peak using the extrema around the R peak on the ECG
lTpeak = length(Tpeak_conv);
for iR = firstiR:length(pkR)
    if (pkR(iR)+round(aftR*Fs)) <= lTpeak
        [~, tmp] = max(Tpeak_conv(pkR(iR)+round(RTlag*Fs) : pkR(iR)+round(aftR*Fs)));
        pkT(iR) = tmp + pkR(iR)+round(RTlag*Fs);
    else
        pkT(iR) = NaN;
    end
end

% COMPUTE CARDIAC TIMINGS AND AMPLITUDE
% =====================================
% convert to beat per minute
RR    = (diff(pkR(:))/Fs)';
PR    = (pkR(:) - pkP(:))/Fs;
QR    = (pkR(:) - pkQ(:))/Fs;
RS    = (pkS(:) - pkR(:))/Fs;
QS    = (pkS(:) - pkR(:))/Fs;
RT    = (pkT(:) - pkR(:))/Fs;

Ramp = zeros(1, length(pkR));
ind = ~any([pkS(:)==0, pkQ(:)==0], 2) ;
Ramp(ind)  = ECG_f(pkR(ind)) - 0.5*(ECG_f(pkS(ind)) + ECG_f(pkQ(ind)));

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function pk = FindECGPeaks(ECG, Fs)
% Find the R peaks in an ECG trace (a vector time series). Return the 
% sample index list for each detected R peaks.
%
% Note: to be accurate, it is better to de-trend the ECG before entering
% this function (e.g. use a highpass filter).
% Method: define a threshold on the ECG trace and look for peak higher than
% this threshold (the R peaks).
%   DEFINITION OF THE THRESHOLD
%       - compute the mean beat rate (by locating the maximal power in the
%       [0.8 - 2.2] frequency band of the ECG power spectrum
%       - use this beat rate to epoch the signal using a constant epoch
%       size of 2*T (T being the above defined beat period)
%           NOTE: one could use arbitrarily T = 1s
%       - the maximum and minimum are searched in each epoch
%       - the mean lower and higher values are defined over epochs and the 
%       threshold set at the highest fourth of this range
%
% Usage: pk = FindECGPeaks(ECG, Fs)
%   ECG: preprocess ECG time series
%    Fs: time series sampling frequency
%    pk: sample indices of R peaks
%
% Florent Meyniel 2011-06-01


% GET THE ECG MEAN HEART BEAT RATE
% ================================

% compute a fourrier transform
L    = length(ECG);
NFFT = 2^nextpow2(L);
Y    = fft(ECG, NFFT)/L;             % abs(Y), i.e. the modulus of Y, is the power for each frequency
f    = Fs/2*linspace(0,1,NFFT/2+1);  % frequency range

% restrict to a specific frequency band (48 - 132 beats per minute)
fmin = 0.8;
fmax = 2.2;
ind  = find(f<fmax & f>fmin);
f_r  = f(ind);
Y_r  = abs(Y(1:NFFT/2+1));
Y_r  = Y_r(ind);

% binning the frequency per 0.05 Hz
binsize = 0.05;
binlist = min(f_r):binsize:max(f_r);
binY    = zeros(1,length(binlist)-1);
for i = 1:length(binlist)-1
    binY(i) = mean(Y_r(f_r>=binlist(i) & f_r<binlist(i+1)));
end 
    
% find peak frequency
[~, ind] = max(binY);
meanHBR  = binlist(ind);

% DEFINE THE MEAN UPPER AND LOWER BOUNDS
% ======================================

% look for extrema in the signal epoched with a period equal to twice the 
% mean HBR, to ensure that there is at least one R peak per epoch.
indlist    = 1:length(ECG);
binsize    = 2/meanHBR*Fs;
binlist    = 1:binsize:length(ECG);
UpperBound = zeros(1,length(binlist)-1);
LowerBound = zeros(1,length(binlist)-1);
for i = 1:length(binlist)-1
    list = indlist>=binlist(i) & indlist<binlist(i+1);
    UpperBound(i) = max(ECG(list));
    LowerBound(i) = min(ECG(list));
end 

% SEARCH FOR R PEACKS
% ===================

% define a search threshold, above 3/4 of the signal range (between its
% bounds)
thd = mean(UpperBound)- (mean(UpperBound) - mean(LowerBound))/3;

% first pass: look for all local peaks above threshold
[~, pk] = findpeaks(ECG, 'MINPEAKHEIGHT', thd);

% Second pass: remove spurious peak
T = diff(pk);
Tmin = mean(1/meanHBR)/2;
pk = pk(~(T<Tmin));

end
