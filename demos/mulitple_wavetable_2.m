%% Multiple wavetable synthesis
% Demonstrates the case of loading wavetable data from audio files.

clear; close all;

% Constants
% Sample rate.
Fs = 44100;
% Output frequency.
F0 = 142.77;
% Output duration.
outDurationS = 4;
% Output amplitude.
outAmp = 1.;
% Wavetable length.
wtLength = 2^9;
% Max output samples to plot.
maxOutPlot = 1000;

% Load some audio into wavetables.
x1 = audioread('./wavetables/vox_wt.wav');
x2 = audioread('./wavetables/violin_wt.wav');
% Resample the wavetables to ensure that they're they contain the same number of
% samples.
wt1 = resample(x1, wtLength, length(x1));
wt2 = resample(x2, wtLength, length(x2));

figure('Name', 'Wavetables', 'Position', [100, 500, 600, 500]), ...
    subplot(211), ...
    plot(1:wtLength, wt1, 'r-'), ...
    title('Wavetable 1'), ...
    xlim([1, wtLength]), ...
    ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('sample index');
    subplot(212), ...
    plot(1:wtLength, wt2, 'r-'), ...
    title('Wavetable 2'), ...
    xlim([1, wtLength]), ...
    ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('sample index');

fig2 = figure('Name', 'Wavetable transition', 'Position', [1000, 500, 600, 500]);

% Create output placeholder.
y = zeros(Fs * outDurationS, 1);

% Calculate the number of samples per period of the wavetable to produce the
% desired frequency.
sampsPerPeriod = Fs / F0;
wtStepsPerSample = wtLength / sampsPerPeriod;

% Placeholder for the transitional wavetable samples.
wt = zeros(2, 1);

% Transition linearly between the wavetables while writing to the output vector.
for n=1:(Fs * outDurationS)    
    % Calculate the fractional sample index in the wavetable.
    wtIndex = mod(wtStepsPerSample * (n - 1), wtLength) + 1;
    % Calculate the magnitudes for samples either side of the fractional index.
    magnitude2 = rem(wtIndex, 1);
    magnitude1 = 1 - magnitude2;
    
    % Wrap the wavetable index as necessary.
    prevIndex = floor(wtIndex);
    nextIndex = ceil(wtIndex);
    if prevIndex == 0
        prevIndex = wtLength;
    end
    if nextIndex > wtLength
        nextIndex = 1;
    end
    
    % Calculate the transitional wavetable samples.
    transition = n / (Fs * outDurationS);
    wt(1) = (1 - transition) * wt1(prevIndex) + transition * wt2(prevIndex);
    wt(2) = (1 - transition) * wt1(nextIndex) + transition * wt2(nextIndex);
    
    % Compose the output sample from summed, weighted transitional samples.
    y(n) = outAmp * ( ...
        magnitude1 * wt(1) + ...
        magnitude2 * wt(2) ...
    );

    % Make a cool plot of the transition (slow for factors of Fs; doesn't work
    % well for large F0).
    if mod(n, sampsPerPeriod) < 0.1
        figure(fig2),
        plot(y(n - floor(sampsPerPeriod) + 1: n), 'k.-'), ...
            ylim([-1, 1]), ...
            title('Wavetable transition'), ...
            drawnow limitrate;
    end
end

sound(y, Fs);
% Plot beginning/end of output against time.
figure( ...
        'Name', 'Output', ...
        'Position', [500 50 750 900] ...
    ),
    subplot(311), ...
    plot(linspace(0, maxOutPlot / Fs, maxOutPlot)', y(1:maxOutPlot)), ...
    title(sprintf('Output waveform (first %d samples)', maxOutPlot)), ...
    ylim([-1, 1]), ...
    ylabel('amp.'), ...
    xlabel('time (ms)'), ...
    subplot(312), ...
    plot(linspace(0, maxOutPlot / Fs, maxOutPlot)', y(length(y) - maxOutPlot + 1:end)), ...
    title(sprintf('Output waveform (last %d samples)', maxOutPlot)), ...
    ylim([-1, 1]), ...
    ylabel('amp.'), ...
    xlabel('time (ms)'), ...
    % Plot spectrogram.
    subplot(313), ...
    spectrogram(y, 512, 64, 512, Fs, 'yaxis'), ...
    ylim([0, 10]);

tfPlot(y, Fs);