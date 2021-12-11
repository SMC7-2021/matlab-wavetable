%% Single wavetable synthesis
% Variable frequency via rounding/tructation.
% Note that, especially for high frequencies, using the sawtooth/square
% wavetable, aliasing artifacts are audible in the output.

clear; close all;

% Constants.
% Mode: 'round' or 'truncate'
mode = 'round';
% Wavetable type: 'sine', 'square', or 'saw'
wtType = 'square';
% Sample rate.
Fs = 44100;
% Output frequency.
F0 = 501;
% Output duration.
outDurationS = 2;
% Output amplitude.
outAmp = 1.;
% Wavetable length.
wtLength = 2^9;
% Max output samples to plot.
maxOutPlot = 500;

% Plot setup.
figure( ...
    'Name', 'Single wavetable', ...
    'Position', [500 50 750 900] ...
);

% Create and plot the wavetable.
switch(wtType)
    case 'sine'
        wt = sin(linspace(0, 2 * pi, wtLength)');
    case 'saw'
        wt = linspace(-1, 1, wtLength);
    case 'square'
        wt = square(linspace(0, 2 * pi, wtLength)');
end

subplot(311), ...
    plot(1:wtLength, wt, 'r.'), ...
    title('Wavetable'), ...
    xlim([1, wtLength]), ...
    ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('sample index');

% Create output placeholder.
y = zeros(Fs * outDurationS, 1);

% Calculate the number of samples per period of the wavetable to produce the
% desired frequency.
sampsPerPeriod = Fs / F0;
wtStepsPerSample = wtLength / sampsPerPeriod;

% Repeatedly copy the wavetable to the output placeholder.
for n=1:(Fs * outDurationS)
    % Calculate the fractional sample index in the wavetable.
    wtIndex = mod(wtStepsPerSample * (n - 1), wtLength) + 1;
    
    % Get the integer sample index as per the mode setting.
    switch(mode)
        case 'round'
            wtIndex = round(wtIndex);
        case 'truncate'
            wtIndex = floor(wtIndex);
    end
    
    % Wrap as necessary.
    if wtIndex == 0
        wtIndex = wtLength;
    end
    if wtIndex > wtLength
        wtIndex = 1;
    end
   
    % Compose the output.
    y(n) = wt(wtIndex);
end

sound(y, Fs);
% Plot beginning of output against time.
subplot(312), ...
    plot(linspace(0, maxOutPlot / Fs, maxOutPlot)', y(1:maxOutPlot)), ...
    title(sprintf('Output waveform (first %d samples)', maxOutPlot)), ...
    ylim([-1, 1]), ...
    ylabel('amp.'), ...
    xlabel('time (ms)');

% Plot spectrogram
subplot(313), ...
    spectrogram(y, 512, 64, 512, Fs, 'yaxis'), ...
    ylim([0, 10]);

% For sine wavetable, inspect the frequency response.
if strcmp(wtType, 'sine')
    figure( ...
        'Name', 'Frequency response', ...
        'Position', [600 50 750 900] ...
    );
    freqz(y);

    % Compare sin wavetable output with 'pure' sine wave:
    figure( ...
        'Name', 'Pure sine', ...
        'Position', [700 50 750 900] ...
    );
    y_ = sin(2 * pi * F0 * linspace(0, outDurationS, Fs * outDurationS));
    % spectrogram( y_, 512, 64, 512, Fs, 'yaxis');
    freqz(y_);
end

tfPlot(y, Fs);
