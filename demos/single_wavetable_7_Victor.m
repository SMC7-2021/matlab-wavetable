%% Single wavetable synthesis
% Variable frequency via sample rate conversion.
% Behaves... weirdly for conversion ratios that require high L,M values. And
% the antialiasing filter seems to be wiped out by downsampling.

clear; close all;

% Constants.
% Wavetable type: 'sine', 'saw', 'square' or 'sample'
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
% Anti-aliasing (low-pass) filter type: 'butter', 'cheby'
aaFilterType = 'butter';
% LPF params
% lpfCutoff = (L * Fs) / (2 * M);
% Cutoff as a proportion of Fs * L...
lpfCutoff = (1/M) * .5;
% [lpfOrder, Wn] = buttord(lpfCutoff * .5, lpfCutoff, 3, 60);
lpfOrder = 6;
% lpfOrder = 5;
% Max output samples to plot.
maxOutPlot = 2^9;

% Set up figure 1.
figure1 = figure( ...
    'Name', 'Single wavetable', ...
    'Position', [400 50 1200 900] ...
);

% Create and plot the wavetable.
switch(wtType)
    case 'sine'
        wt = sin(linspace(0, 2 * pi, wtLength)');
    case 'saw'
        wt = linspace(-1, 1, wtLength)';
    case 'square'
        wt = square(linspace(0, 2 * pi, wtLength)');
    case 'sample'
        x = audioread('./wavetables/vox_wt.wav');
        wt = resample(x, wtLength, length(x));
end

switch(aaFilterType)
    case 'butter'
%         [b, a] = butter(lpfOrder, lpfCutoff/(Fs/2), 'low');
        [b, a] = butter(lpfOrder, lpfCutoff, 'low');
    case 'cheby'
        % 'peak to peak passband ripple', dB
        ripple = 3;
%         [b, a] = cheby1(lpfOrder, ripple, lpfCutoff/(Fs/2), 'low');
        [b, a] = cheby1(lpfOrder, ripple, lpfCutoff, 'low');
end

% Plot the original wavetable
f = linspace(-Fs/2 + 1, Fs/2 - 1, Fs)';
wt1s = repmat(wt, ceil(Fs/wtLength), 1);
wt1s = wt1s(1:Fs);
subplot(521), ...
    plot(1:wtLength, wt, 'r.'), ...
    title('Original wavetable'), ...
    xlim([1, wtLength]), ...
    ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('sample index');
subplot(522), ...
    plot(f, abs(fftshift(fft(wt1s)))), ...
    title('Spectrum of original wavetable');

% Upsample the wavetable
wtUp = L * upsample(wt, L);
FsUp = Fs*L;
fUp = linspace(-FsUp/2 + 1, FsUp/2 - 1, FsUp)';
wtUp1s = repmat(wtUp, ceil(FsUp/length(wtUp)), 1);
wtUp1s = wtUp1s(1:FsUp);
subplot(523), ...
    plot(1:length(wtUp), wtUp), ...
    title('Upsampled wavetable'), ...
%     ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('sample index');
subplot(524), ...
    plot(fUp, abs(fftshift(fft(wtUp1s)))), ...
    title('Spectrum of upsampled wavetable');

% Filter the upsampled wavetable
z = zeros(1, lpfOrder);
wtUpFlt = wtUp;
for i=1:3
    [wtUpFlt, z] = filter(b, a, wtUpFlt, z);
end
wtUpFlt1s = repmat(wtUpFlt, ceil(FsUp/length(wtUpFlt)), 1);
wtUpFlt1s = wtUpFlt1s(1:FsUp);
subplot(525), ...
    plot(1:length(wtUpFlt), wtUpFlt), ...
    title('Filtered upsampled wavetable'), ...
%     ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('sample index');
subplot(526), ...
    plot(fUp, abs(fftshift(fft(wtUpFlt1s)))), ...
    title('Spectrum of filtered upsampled wavetable');

% Downsample the filtered wavetable.
wtDn = downsample(wtUpFlt, M);
FsDn = (Fs*L)/M;
fDn = linspace(-FsDn/2 + 1, FsDn/2 - 1, FsDn)';
wtDn1s = repmat(wtDn, ceil(FsDn/length(wtDn)), 1);
wtDn1s = wtDn1s(1:FsDn);
subplot(527), ...
    plot(1:length(wtDn), wtDn), ...
    title('Downsampled wavetable'), ...
%     ylim([-1.1, 1.1]), ...
    ylabel('amp.'), ...
    xlabel('sample index');
subplot(528), ...
    plot(fDn, abs(fftshift(fft(wtDn1s)))), ...
    title('Spectrum of downsampled wavetable');

wtFinal = wtDn;

% Create output placeholder.
y = zeros(Fs * outDurationS, 1);
% Initialize wavetable read index.
wtIndex = 0;

% Repeatedly copy the wavetable to the output placeholder.
for n=1:(Fs * outDurationS)
    wtIndex = mod(wtIndex, length(wtFinal)) + 1; 
    y(n) = wtFinal(wtIndex);
end

soundsc(y, Fs);
% Plot beginning of output against time.
figure(figure1);
subplot(529), ...
    plot(linspace(0, maxOutPlot / Fs, maxOutPlot)', y(1:maxOutPlot)), ...
    title(sprintf('Output waveform (first %d samples)', maxOutPlot)), ...
%     ylim([-1.25, 1.25]), ...
    ylabel('amp.'), ...
    xlabel('time (ms)');

% Plot spectrogram
subplot(5,2,10), ...
    spectrogram(y, 512, 64, 512, Fs, 'yaxis'), ...
    title('Spectrogram of output'), ...
    ylim([0, 22.05]);

tfPlot(y, Fs, .01);