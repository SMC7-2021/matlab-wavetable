clear; close all; clc;

Fs = 44100;
duration = 1;
F0 = 440;

tic
y = wavetable(Fs, duration, F0, ...
    'WavetableType', 'square', ... % Only 'square' available for now...
    'InterpolationType', 'cubic', ...
    'OutputAmplitude', 1, ...
    'Oversample', 2 ...
);
toc

tfPlot(y, Fs);
soundsc(y, Fs);

maxOutPlot = 1000;
figure( ...
        'Name', 'Output', ...
        'Position', [500 50 750 900] ...
    ),
    subplot(311), ...
    plot(linspace(0, maxOutPlot / Fs, maxOutPlot)', y(1:maxOutPlot)), ...
    title(sprintf('Output waveform (first %d samples)', maxOutPlot)), ...
    ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('time (ms)'), ...
    subplot(312), ...
    plot(linspace(0, maxOutPlot / Fs, maxOutPlot)', y(length(y) - maxOutPlot + 1:end)), ...
    title(sprintf('Output waveform (last %d samples)', maxOutPlot)), ...
    ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('time (ms)'), ...
    subplot(313), ...
    spectrogram(y, 512, 64, 512, Fs, 'yaxis'), ...
    ylim([0, 22]);