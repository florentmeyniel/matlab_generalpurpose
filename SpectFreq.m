function [binY, binlist] = SpectFreq(data, Fs, minF, maxF, binSize, doplot)
% Plot the power spectrum of a vector time series using fast- fourrier
% transform.
%
% Usage: [binY, binlist] = SpectFreq(data, Fs, minF, maxF, binSize, doplot)
%      data: vector time series
%        Fs: data sampling rate
%      minF: minimum frequency to display
%      maxF: maximum frequency to display
%   binSize: bin size (Hz) to bin data (optional)
%    doplot: 1 to plot (default), 0 otherwise
%
% Florent Meyniel 2012-01-06
 

% compute a fourrier transform
L    = length(data);
NFFT = 2^nextpow2(L);
Y    = fft(data, NFFT)/L; % abs(Y), i.e. the modulus of Y is the power for each frequency
f    = Fs/2*linspace(0,1,NFFT/2+1); % frequency range

% check input
if minF>maxF
    error('the specified minimum %d > than maximum %d !!!', minF, maxF)
end
if minF>Fs/2
    error('the specified minimum %d is higher than the Nyquisdt frequency %d', minF, Fs/2)
end
if maxF>Fs/2    
    warning('the specified maximum %d is reset to the Nyquisdt frequency %d', maxF, Fs/2)
    maxF = Fs/2;
end

if exist('binSize', 'var')
    if isempty(binSize)
        doBin = 0;
    else
        if binSize < unique(diff(f))
            warning('no effective binning (bin size: %d; resolution: %d', binSize, unique(diff(f)))
            doBin = 0;
        else
            doBin = 1;
        end
    end
else
    doBin = 0;
end

if ~exist('doplot', 'var')
    doplot = 1;
end

% restrict to a specific frequency band
ind  = find(f<maxF & f>minF);
f_r  = f(ind);
Y_r  = abs(Y(1:NFFT/2+1));
Y_r  = Y_r(ind);

if doBin
    % binning the frequency per 0.1 Hz
    binlist = min(f_r):binSize:max(f_r);
    binY    = zeros(1,length(binlist)-1);
    for i = 1:length(binlist)-1
        binY(i) = mean(Y_r(f_r>=binlist(i) & f_r<binlist(i+1)));
    end
else
    binY    = Y_r;
    binlist = [f_r 0];
end

binlist = binlist(1:end-1);

if doplot
    iFig = figure;
    set(iFig, 'Name', 'Spectral Analysis')
    set(iFig, 'Color', [1 1 1])
    plot(binlist, binY, 'LineWidth', 2)
    ylabel('Power')
    xlabel('Hz')
end

end