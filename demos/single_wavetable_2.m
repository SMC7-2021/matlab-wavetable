%% Single wavetable
% Variable frequency via interpolation.
% Note that, especially for high frequencies, rounding discrepancies introduce
% interpolation noise into the output signal.
% One possible workaround could be to up/down-sample the wavetable to fit the 
% desired frequency before the loop... but this wouldn't be practical in real-
% time.

clear; close all;

% Constants.
% Wavetable type: 'sine' or 'saw'
wtType = 'sine';
% Sample rate.
Fs = 44100;
% Output frequency.
F0 = 600;
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
    
    % Compose the output sample from summed weighted samples.
    y(n) = outAmp * ( ...
        magnitude1 * wt(prevIndex) + ...
        magnitude2 * wt(nextIndex) ...
    );
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
