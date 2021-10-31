clear; close all;

%% Single sawtooth wavetable
% Output frequency can only be varied by changing:
% - the length of the wavetable;
% - the sample rate;

% Constants.
% Sample rate.
Fs = 44100;
% Output duration.
outDurationS = 2;
% Output amplitude.
outAmp = .75;
% Wavetable length.
wtLength = 75;
% Max output samples to plot.
maxOutPlot = 2000;

% Plot setup.
figure( ...
    'Name', 'Basic single wavetable', ...
    'Position', [500 50 750 900] ...
);

% Create and plot the wavetable.
wt = linspace(-1, 1, wtLength);

subplot(311), ...
    plot(1:wtLength, wt, 'r.'), ...
    title('Wavetable'), ...
    xlim([1, wtLength]), ...
    ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('sample index');

% Create output placeholder.
y = zeros(Fs * outDurationS, 1);

% Repeatedly copy the wavetable to the output placeholder.
for n=1:(Fs * outDurationS)
    wtIndex = mod(n-1, wtLength) + 1;
    y(n) = outAmp * wt(wtIndex);
end

sound(y, Fs);
% Plot beginning of output against time.
subplot(312), ...
    plot(linspace(0, maxOutPlot / Fs, maxOutPlot)', y(1:maxOutPlot)), ...
    title(sprintf('Output waveform (first %d samples)', maxOutPlot)), ...
    ylim([-1, 1]), ...
    ylabel('amp.'), ...
    xlabel('time (ms)');

subplot(313), ...
    spectrogram(y, 512, 64, 512, Fs, 'yaxis');
    