clc; clear; close all; 
addpath('../helpers')

% Sampling rate.
Fs = 44100;
% Output frequency.
F0 = 440;

% F0 = linspace(66, 2027, Fs * outDurationS)';
% Output duration.
outDurationS = 0.01;
% Wavetable length.
wtLength = 2^7;
% Target samlping rate, required to reproduce requested F0.
targetFs = F0 * wtLength;
targetWtLength = Fs / F0;

wt = sin(linspace(0, (2*pi) - (2*pi)/wtLength, wtLength)');
% wt = square(linspace(0, (2*pi) - (2*pi)/wtLength, wtLength)');
% wt = square(linspace(0, (2*pi), wtLength)');

% Calculate the number of samples per period of the wavetable to produce the
% current frequency.
sampsPerPeriod = Fs / F0;
wtStepsPerSample = 1.2;

%% Drop-sample interpolation

y = zeros(wtLength, 1);
readindex = 1;
for n=1:wtLength
    
    y(n) = wt(readindex);
    readindex = mod(readindex + wtStepsPerSample, wtLength);
    readindex = floor(readindex);
    if readindex == 0
        readindex = length(wt);
    end

end
%tfPlot(y, Fs, .005);
%figure; snr(y, Fs);
%figure; spectrogram(y, 512, 64, 512, Fs, 'yaxis');

%% linear interpolation 
%wtStepsPerSample = wtLength / sampsPerPeriod;
readindex = 1;
y2 = zeros(wtLength, 1);

for n=1:wtLength

    alphaNext = mod(readindex, 1);
    alphaPrev = 1 - alphaNext;

    prevIndex = floor(readindex);
    nextIndex = prevIndex + 1;
    
    if prevIndex == 0
        prevIndex = length(wt);
    end
    if nextIndex == length(wt)
        nextIndex = 1;
    end

    y2(n) = alphaPrev * wt(prevIndex) + alphaNext * wt(nextIndex);

    readindex = mod(readindex + wtStepsPerSample, wtLength);
    
end

%tfPlot(y2, Fs, .005);
%figure; snr(y2, Fs);
%figure; spectrogram(y2, 512, 64, 512, Fs, 'yaxis');

%% Cubic Lagrange
%wtStepsPerSample = wtLength / sampsPerPeriod;
readindex = 1;
y3 = zeros(wtLength,1);

for n=1:wtLength
    % Get the indices of the samples we need from the wavetable.
    % Start by getting the index of the ith sample.
    i = floor(readindex);
    % Create a vector of the four sample indices we need for cubic interpolation.
    readIndices = [i-1, i, i+1, i+2];
    % Wrap the indices as required.
    for i=1:length(readIndices)
        if readIndices(i) < 1
            readIndices(i) = length(wt) + readIndices(i);
        elseif readIndices(i) > length(wt)
            readIndices(i) = 1 + mod(readIndices(i), length(wt));
        end
    end
    
    % Get the fractional part of the wt read index.
    alpha = mod(readindex, 1);
    % Calculate the sample amplitude coefficients.
    % coeffs(1) is the coefficient for the (n-1)th sample
    % coeffs(2) => nth sample
    % coeffs(3) => (n+1)th sample
    % coeffs(4) => (n+2)th sample
    coeffs = [ ...
        -alpha * (alpha - 1) * (alpha - 2) / 6, ...
        (alpha - 1) * (alpha + 1) * (alpha - 2) / 2, ...
        -alpha * (alpha + 1) * (alpha - 2) / 2, ...
        alpha * (alpha + 1) * (alpha - 1)/6 ...
    ];

%     y3(n) = coeffs(1)*wt(readIndices(1)) + ...
%         coeffs(2)*wt(readIndices(2)) + ...
%         coeffs(3)*wt(readIndices(3)) + ...
%         coeffs(4)*wt(readIndices(4));
    % One-liner equivalent to the above:
    y3(n) = coeffs * wt(readIndices);

    readindex = mod(readindex + wtStepsPerSample, wtLength);
    
end
% snr(sin(2*pi*440*linspace(0,outDurationS,Fs*outDurationS)));
%tfPlot(y3, Fs, .005);
%figure; snr(y3, Fs);
%figure; spectrogram(y3, 512, 64, 512, Fs, 'yaxis');

%%
figure
hold on
stem(wt)
plot(y)
plot(y2)
plot(y3)
legend('Original Wavetable', 'Drop-Zero','Linear','Cubic')
grid on
set(gcf, 'color', 'w');    