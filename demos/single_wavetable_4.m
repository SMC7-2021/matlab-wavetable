%% Single wavetable synthesis
% Variable frequency via interpolation.
% Attempted anti-aliasing filter.
% Doesn't sound pretty, especially for high F0.

clear; close all;

% Constants.
% Wavetable type: 'sine', 'saw', or 'square'
wtType = 'square';
% Sample rate.
Fs = 44100;
% Output frequency.
F0 = 2532;
% Output duration.
outDurationS = 2;
% Output buffer size.
bufferSize = 2^10;
% Output amplitude.
outAmp = 1.;
% Wavetable length.
wtLength = 2^11;
% Max output samples to plot.
maxOutPlot = 2^9;
% Anti-aliasing (low-pass) filter type: 'butter', 'cheby'
aaFilterType = 'cheby';
% LPF params
lpfCutoff = Fs * .4;
lpfOrder = 8;

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

switch(aaFilterType)
    case 'butter'
        [b, a] = butter(lpfOrder, lpfCutoff/(Fs/2), 'low');
    case 'cheby'
        % 'peak to peak passband ripple'
        ripple = 3;
        [b,a] = cheby1(lpfOrder, ripple, lpfCutoff/(Fs/2), 'low');
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
buffer = zeros(bufferSize, 1);
bufferIndex = 1;

% Calculate the number of samples per period of the wavetable to produce the
% desired frequency.
sampsPerPeriod = Fs / F0;
wtStepsPerSample = wtLength / sampsPerPeriod;

% Placeholder for final filter conditions.
z = zeros(lpfOrder, 1);

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
    if prevIndex == 0 || prevIndex > wtLength
        prevIndex = wtLength;
    end
    if nextIndex > wtLength
        nextIndex = 1;
    end
    
    % Compose the output sample from summed weighted samples.
    buffer(bufferIndex) = outAmp * ( ...
        magnitude1 * wt(prevIndex) + ...
        magnitude2 * wt(nextIndex) ...
    );

    bufferIndex = bufferIndex + 1;
    if bufferIndex > bufferSize
        % End of buffer reached: low-pass the buffer and write it to output.
        [filtered, z] = filter(b, a, buffer, z);
        y(n-bufferSize+1:n) = filtered;
        bufferIndex = 1;
    end
end

sound(y, Fs);
% Plot beginning of output against time.
subplot(312), ...
    plot(linspace(0, maxOutPlot / Fs, maxOutPlot)', y(1:maxOutPlot)), ...
    title(sprintf('Output waveform (first %d samples)', maxOutPlot)), ...
    ylim([-1.25, 1.25]), ...
    ylabel('amp.'), ...
    xlabel('time (ms)');

% Plot spectrogram
subplot(313), ...
    spectrogram(y, 512, 64, 512, Fs, 'yaxis'), ...
    ylim([0, 22]);

% % For sine wavetable, inspect the frequency response.
% if strcmp(wtType, 'sine')
%     figure( ...
%         'Name', 'Frequency response', ...
%         'Position', [600 50 750 900] ...
%     );
%     freqz(y);
% 
%     % Compare sin wavetable output with 'pure' sine wave:
%     figure( ...
%         'Name', 'Pure sine', ...
%         'Position', [700 50 750 900] ...
%     );
%     y_ = sin(2 * pi * F0 * linspace(0, outDurationS, Fs * outDurationS));
%     % spectrogram( y_, 512, 64, 512, Fs, 'yaxis');
%     freqz(y_);
% end

tfPlot(y, Fs);
