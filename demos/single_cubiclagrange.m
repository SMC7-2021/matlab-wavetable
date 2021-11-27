% Constants.

clc; clear all; 

%%

% Wavetable type: 'sine', 'saw', or 'square'
wtType = 'square';
% Sampling rate.
Fs = 44100;
% Output frequency.
F0 = 440;

% F0 = linspace(66, 2027, Fs * outDurationS)';
% Output duration.
outDurationS = 0.04;
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

wtStepsPerSample = wtLength / sampsPerPeriod;
readindex = 1;
y2 = zeros(Fs * outDurationS,1);

for n=1:(Fs * outDurationS)-1

    readindex = mod(readindex + wtStepsPerSample,wtLength) +1 ;

    y(n) = wt(floor(readindex));




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
wtStepsPerSample = wtLength / sampsPerPeriod;
readindex = 1;
y3 = zeros(Fs * outDurationS,1);

for n=1:Fs * outDurationS

    wtIndex = mod(wtStepsPerSample * (n - 1), wtLength) + 1;

    alpha = rem(wtIndex, 1);

    test(n)=alpha;
    
    if (alpha ~= 0)

        l0= (-alpha*(alpha-1)*(alpha-2)) /6;
        l1= (alpha-1)*(alpha+1)*(alpha-2) /2;
        l2= (-alpha*(alpha+1)*(alpha-2)) /2;
        l3= (alpha*(alpha+1)*(alpha-1)) /6;
    end
    
    readindex = floor(wtIndex);
    x0 = mod(readindex + wtLength - 1,wtLength); % Length
    x1 = readindex;
    x2 = mod(readindex +1, wtLength+1);
    x3 = mod(readindex +2, wtLength+1);
    if (x0 == 0) 
        x0 = wtLength; 
    end
        if (x3 == 0) 
        x3 = 1; 
    end
    if (x2 == 0) 
        x2 = 1; x3 = 2;
    end

    if (alpha ~= 0)
    y3(n) = ...
        l0*wt(x0) + ...
        l1*wt(x1) + ...
        l2*wt(x2) + ...
        l3*wt(x3);
    else 
        y3(n) = wt(readindex);
    end
   % y3 = y3';
end

%%

%tfPlot(y,Fs,.01)
%tfPlot(y2,Fs,.01)
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

