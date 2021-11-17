%% Single wavetable synthesis
% Variable frequency via downsampling.
% Attempted anti-aliasing filter.

clear; close all;

% Constants.
% Wavetable type: 'sine', 'saw', or 'square'
wtType = 'square';
% Sampling rate.
Fs = 44100;
% Output frequency.
F0 = 440;
% Output duration.
outDurationS = 2;
% Wavetable length.
wtLength = 2^11;
% Target samlping rate, required to reproduce requested F0.
targetFs = F0 * wtLength;
targetWtLength = Fs / F0;
% Upsampling/downsampling factors
[L, M] = getResamplingFactors(wtLength, targetWtLength);
M = round(wtLength / (Fs / F0));
% M = 9; L = 3;
% Anti-aliasing (low-pass) filter type: 'butter', 'cheby'
aaFilterType = 'cheby';
% LPF params
lpfCutoff = Fs * .35; % Fs * M * .35;
lpfOrder = 10;
% Output amplitude.
outAmp = 1.;
% Max output samples to plot.
maxOutPlot = 2^9;

% Plot setup.
figure1 = figure( ...
    'Name', 'Single wavetable', ...
    'Position', [500 50 750 900] ...
);

% Create and plot the wavetable.
switch(wtType)
    case 'sine'
        wt = sin(linspace(0, 2 * pi, wtLength)');
    case 'saw'
        wt = linspace(-1, 1, wtLength)';
    case 'square'
        wt = square(linspace(0, 2 * pi, wtLength)');
end

switch(aaFilterType)
    case 'butter'
%         [b, a] = butter(lpfOrder, lpfCutoff/((Fs*M)/2), 'low');
%         [b, a] = butter(lpfOrder, min([1/M, 1/L]), 'low');
        [b, a] = butter(lpfOrder, lpfCutoff/(Fs/2), 'low');
    case 'cheby'
        % 'peak to peak passband ripple', dB
        ripple = 3;
%         [b, a] = cheby1(lpfOrder, ripple, min([1/M, 1/L]), 'low');
%         [b, a] = cheby1(lpfOrder, ripple, .05, 'low');
        [b, a] = cheby1(lpfOrder, ripple, lpfCutoff/(Fs/2), 'low');

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

% Resample the wavetable
wtUpsampled = wt; %upsample(wt, L);
wtRepeated = repmat(wtUpsampled, 3, 1);
wtRepeatedDownsampled = downsample(wtRepeated, M);
wtFiltered = filter(b, a, wtRepeatedDownsampled);
wtFinal = wtFiltered(end - length(wtRepeatedDownsampled) / 3 + 1: end);
wtDownsampled = wtFinal;

% z = zeros(lpfOrder, 1);
% wtDownsampled = downsample(wt, M);
% wtFiltered = wtDownsampled;
% for i=1:2
%     [wtFiltered, z] = filter(b, a, wtFiltered, z);
% end
% wtFinal = wtFiltered;

% Final filter conditions
% z = zeros(lpfOrder, 1);
% [wtFiltered, z] = filter(b, a, wtRepeated, z);
% [wtFiltered, z] = filter(b, a, wtUpsampled, z);
% wtDownsampled = downsample(wtFiltered(end-length(wtUpsampled)+1:end), M);
% wtDownsampled = downsample(wtFiltered, M);
finalWtLength = length(wtFinal);

figure2 = figure( ...
    'Name', 'Single wavetable', ...
    'Position', [1000 50 750 900] ...
);

subplot(411), ...
    plot(wtDownsampled), ...
    title('Downsampled');
subplot(412), ...
    plot(wtFiltered), ...
    title('Filtered');
subplot(413), ...
    plot(wtFinal), ...
    title('Final');
subplot(414), ...
    spectrogram(wtFiltered, 64, 8, 128, Fs, 'yaxis');

wtIndex = 0;


% Repeatedly copy the wavetable to the output placeholder.
for n=1:(Fs * outDurationS)
    wtIndex = mod(wtIndex, finalWtLength) + 1;
    y(n) = wtFinal(wtIndex);
end

soundsc(y, Fs);
% Plot beginning of output against time.
figure(figure1);
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
