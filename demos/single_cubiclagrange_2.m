% Constants.

clc; clear all; 

% Wavetable type: 'sine', 'saw', or 'square'
wtType = 'square';
% Sampling rate.
Fs = 44100;
% Output duration.
outDurationS = 9;
% Output frequency.
%F0 = 440;

 F0 = linspace(66, 7902, Fs * outDurationS)';

% Wavetable length.
wtLength = 2^9;
% Target samlping rate, required to reproduce requested F0.
targetFs = F0 * wtLength;
targetWtLength = Fs / F0;

wt = square(linspace(0, 2 * pi, wtLength)');

readindex = 1;



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

    wtIndex = mod(wtStepsPerSample * (n - 1), wtLength) + 1;

    alphaNext = rem(wtIndex, 1);
    alphaPrev = 1 - alphaNext;

    prevIndex = mod(floor(readindex), 1);
    nextIndex = prevIndex + 1;

        % Wrap the wavetable index as necessary.
    prevIndex = floor(wtIndex);
    nextIndex = ceil(wtIndex);

    if prevIndex == 0 || prevIndex > wtLength
        prevIndex = wtLength;
    end
    if nextIndex > wtLength
        nextIndex = nextIndex-wtLength;
    end
    
    y2(n) = alphaPrev * wt(prevIndex) + alphaNext * wt(nextIndex);

end

%% cubic lagrange interpolation 

readindex = 1;
y2 = zeros(Fs * outDurationS,1);

for n=1:Fs * outDurationS

    % Calculate the number of samples per period of the wavetable to produce the
    % current frequency.
    sampsPerPeriod = Fs / F0(n);
    wtStepsPerSample = wtLength / sampsPerPeriod;

   % wtStepsPerSample = floor(wtStepsPerSample);

    wtIndex = mod(wtStepsPerSample * (n - 1), wtLength) + 1;

    alpha = rem(wtIndex, 1);
    
    if (alpha ~= 0)

        l0= (-alpha*(alpha-1)*(alpha-2)) /6;
        l1= (alpha-1)*(alpha+1)*(alpha-2) /2;
        l2= (-alpha*(alpha+1)*(alpha-2)) /2;
        l3= (alpha*(alpha+1)*(alpha-1)) /6;
    end
    
    readindex = floor(wtIndex);
    x0 = readindex-1;
    x1 = readindex;
    x2 = readindex+1;
    x3 = readindex+2;
    if (x0 < wtLength) x0 = wtLength; end
    if (x2 > wtLength) x2 = x2 - wtLength; end
    if (x3 > wtLength) x3 = x3 - wtLength; end

    if (alpha ~= 0)
    y3(n) = ...
        l0*wt(x0) + ...
        l1*wt(x1) + ...
        l2*wt(x2) + ...
        l3*wt(x3);
    else 
        y3(n) = wt(readindex);
    end
end

%%

tfPlot(y,Fs,.01)
tfPlot(y2,Fs,.01)
tfPlot(y3,Fs,.01)

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

