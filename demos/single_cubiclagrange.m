% Constants.

clc; clear all; 

% Wavetable type: 'sine', 'saw', or 'square'
wtType = 'square';
% Sampling rate.
Fs = 44100;
% Output frequency.
F0 = 440;

% F0 = linspace(66, 2027, Fs * outDurationS)';
% Output duration.
outDurationS = 0.5;
% Wavetable length.
wtLength = 2^9;
% Target samlping rate, required to reproduce requested F0.
targetFs = F0 * wtLength;
targetWtLength = Fs / F0;

wt = sin(linspace(0, 2 * pi, wtLength)');

readindex = 1;

% Calculate the number of samples per period of the wavetable to produce the
% current frequency.
sampsPerPeriod = Fs / F0;
wtStepsPerSample = wtLength / sampsPerPeriod;

wtStepsPerSample = floor(wtStepsPerSample);

y = zeros((Fs * outDurationS)-1,1);

%% Drop-sample interpolation
for n=1:(Fs * outDurationS)-1


    y(n) = wt(readindex);

    readindex = mod(readindex + wtStepsPerSample,wtLength) +1 ;


end

%% linear interpolation 
wtStepsPerSample = wtLength / sampsPerPeriod;
readindex = 1;
y2 = zeros(Fs * outDurationS,1);

for n=1:Fs * outDurationS

    alphaNext = mod(readindex, 1);
    alphaPrev = 1 - alphaNext;

    prevIndex = mod(floor(readindex),;
    nextIndex = prevIndex + 1;

    y2(n) = alphaPrev * wt(prevIndex) + alphaNext * wt(nextIndex);

    readindex = mod(readindex + wtStepsPerSample, wtLength);
    
end

%tfPlot(y,Fs,.1)

%%
function tfPlot(x, Fs, varargin)
%TFPLOT Create and display a time/frequency plot
    N = length(x);
    t = (1:N)/Fs;
    
    figure('Position', [500, 300, 1000, 400]);
    
    subplot(121),...
        plot(t, x),...
        grid on,...
        xlabel('Time (sec)');
    
    % Limit the x-range for the time plot if a limit (in seconds) has been 
    % specified.
    if nargin == 3 && isnumeric(varargin{1})
        xlim([0, varargin{1}]);
    end
    
    w = chebwin(length(x), 200);
    X = abs(fft(w.*x, Fs));
    X = db(X/max(X));
    subplot(122),...
%         semilogx(linspace(0, Fs-1, Fs), X), ...
        plot(linspace(0, Fs-1, Fs), X, 'LineWidth', 1), ...
        axis([0 Fs/2 -100 5]), ...
        xlabel('Freq (kHz)'),...
        ylabel('Magnitude (dB)'), ...
        title('Spectrum'),...
        set(gca,'XTick',[0 5e3 10e3 15e3 20e3]), ...
        set(gca,'XTickLabel',[{'0'}, '5', '10' '15' '20']), ...
        set(gca,'YTick',-100:20:0), ...
        grid on;
end

