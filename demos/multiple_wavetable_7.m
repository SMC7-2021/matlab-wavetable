%% Precalculated wavetables

% - Load wavetables
% - For each wavetable, create frequency-specific mipmaps
%   - Specify frequency ranges (8ve, 1/3 8ve...)
%   - For each frequency range
%       - Calculate highest allowed frequency component (1)
%       - Compute FFT of wavetable
%       - Attenuate frequency content above limit calculated at (1)
%       - compute IFFT
% - Generate output
%   - For F0, select the appropriate mipmap to use for output.

clear; close all; clc;

% Load some audio samples to use as wavetables.
x1 = audioread('./wavetables/vox_wt.wav');
x2 = audioread('./wavetables/violin_wt.wav');
% Sample rate.
Fs = 44100;
% Output duration.
outDurationS = 1;
% Output amplitude.
outAmp = 1.;
% Wavetable length.
wtLength = 2^11;
shuffleWavetables = false;
% Sample interpolation type: 'truncate', 'linear', 'cubic'
interpolationType = 'cubic';
% Frequency output type: 'sweep', 'sweepmidi', 'fixed440', 'fixed1k', 'fixed2k'
outputType = 'fixed2k';
% Max output samples to plot.
maxOutPlot = 1000;
% (Set up a rate-limiter for the wavetable transition plot.)
transitionPlotIndex = 0;
transitionPlotRateLimit = 200;
doTransitionPlot = true;
doTfPlot = false;

% Create a vector of output frequencies.
switch outputType
    case 'fixed1k'
        F0 = linspace(1000.1, 1000.1, Fs * outDurationS)';
        doTfPlot = true;
    case 'fixed2k'
        F0 = linspace(2000, 2000, Fs * outDurationS)';
        doTfPlot = true;
    case 'sweep'
        outDurationS = 20;
        F0 = linspace(20, 10000, Fs * outDurationS)';
    case 'sweepmidi'
        outDurationS = 20;
        F0 = linspace(midi2hz(1), midi2hz(127), Fs * outDurationS)'; 
    otherwise % fixed440
        F0 = linspace(440, 440, Fs * outDurationS)';
        doTfPlot = true;
end

% Create an array of source wavetables. (Uncomment a statement to enable a
% wavetables; uncomment multiple statements to enable morphing.)
sourceWavetables = {
    % A sampled wavetable.
%     resample(x1, wtLength, length(x1));
    % A sine wavetable.
%     sin(linspace(0, (2*pi) - (2*pi)/wtLength, wtLength)');
    % A discontinuous sine wavetable.
%     sin(linspace(0, 2.2*pi, wtLength)');
    % A square wavetable.
    square(linspace(0, (2*pi) - (2*pi)/wtLength, wtLength)');
    % Another sampled wavetable.
%     resample(x2, wtLength, length(x2));
    % A sawtooth wavetable.
%     sawtooth(linspace(0, (2*pi) - (2*pi)/wtLength, wtLength)');
};
if shuffleWavetables
    sourceWavetables = sourceWavetables(randperm(length(sourceWavetables)));
end

% Plot the source wavetables.
figure('Name', 'Wavetables', 'Position', [100, 500, 600, 900]);

for i=1:length(sourceWavetables)
    subplot(length(sourceWavetables), 1, i), ...
        plot(1:wtLength, sourceWavetables{i}, 'r-'), ...
        title(sprintf('Wavetable %d', i)), ...
        xlim([1, wtLength]), ...
        ylim([-1.1, 1.1]), ...
        ylabel('amp.'), ...
        xlabel('sample index');
end

% Create a placeholder for the wavetable mipmaps.
wavetables = cell(1, length(sourceWavetables));

% Create a vector of frequencies to serve as the basis for the wavetable mipmaps.
% Start from the 'fundamental' frequency of a wavetable of length 'wtLength' at
% sampling rate Fs (44100/2048 = 21.5332…) and proceed in octaves.
% NB, octave spacing leaves audible gaps in the output spectrum, especially when
% sweeping. 1/3 octave might work better.
basisF0s = (Fs/wtLength) * 2.^(0:9);

% Set up wavetable mipmaps.
for i=1:length(sourceWavetables)
    wt = sourceWavetables{i};
    x = zeros(wtLength, length(basisF0s));
    
    % For each fundamental frequency, compute a mipmap of the wavetable 
    % bandlimited to Nyqvist.
    for f = 1:length(basisF0s)
        x(:, f) = computeMipmap(wt, Fs, basisF0s(f));
    end
    
    wavetables{i} = x;
end

% Placeholder for the wavetable morph samples.
wt = zeros(2, 1);
% Initialize the wavetable sample index.
wtSampIndex = 1;
% Create output placeholder.
y = zeros(Fs * outDurationS, 1);

if doTransitionPlot
    fig2 = figure('Name', 'Wavetable transition', 'Position', [1000, 500, 600, 500]);
end

% Generate output.
for n=1:Fs*outDurationS
    % Calculate how far through the output we are.
    outputDelta = n / (Fs * outDurationS);
    % Multiply this by the number of regions between wavetables (e.g. 5 
    % wavetables, 4 regions) to work out how far we are between the current 
    % wavetable pair.
    wtDelta = outputDelta * (length(wavetables) - 1);
    % The integer part is the index of the previous wavetable.
    wtIndex = floor(wtDelta) + 1;
    % The fractional part is the morph amount between the pair of wavetables.
    morph = rem(wtDelta, 1);
    % Calculate the index of the next wavetable.
    nextWtIndex = wtIndex + 1;
    % At the very end of output, wtDelta == length(wavetables) so nextWtIndex 
    % will exceed the length of the wavetable array; so check for this.
    if nextWtIndex > length(wavetables)
        nextWtIndex = length(wavetables) - 1;
    end
    if nextWtIndex == 0
        nextWtIndex = 1;
    end
    
    % Determine, based on the current output frequency, which mipmap to use.
    mipmap = length(basisF0s);
    for f=length(basisF0s):-1:1
        if F0(n) > basisF0s(f)
            mipmap = f;
            break
        end
    end

    % Calculate the transitional wavetable samples.
    switch interpolationType
        case 'truncate'
            wt(1) = interpolateZeroth(wavetables{1, wtIndex}(:, mipmap), wtSampIndex);
            wt(2) = interpolateZeroth(wavetables{1, nextWtIndex}(:, mipmap), wtSampIndex);
        case 'linear'
            wt(1) = interpolateLinear(wavetables{1, wtIndex}(:, mipmap), wtSampIndex);
            wt(2) = interpolateLinear(wavetables{1, nextWtIndex}(:, mipmap), wtSampIndex);
        otherwise % cubic
            wt(1) = interpolateCubic(wavetables{1, wtIndex}(:, mipmap), wtSampIndex);
            wt(2) = interpolateCubic(wavetables{1, nextWtIndex}(:, mipmap), wtSampIndex);
    end

    % Calculate the output sample.
    y(n) = outAmp * (...
        (1 - morph) * wt(1) + ...
        morph * wt(2) ...
    );

    % Calculate the number of samples per period of the wavetable to produce the
    % current frequency.
    sampsPerPeriod = Fs / F0(n);
    wtStepsPerSample = wtLength / sampsPerPeriod;
    
    % (Previous index just used for plotting.)
    prevWtSampIndex = wtSampIndex;
    
    % Upadate the wavetable sample index for the next iteration.
    wtSampIndex = mod(wtSampIndex + wtStepsPerSample, wtLength);

    % Make a cool plot of the wavetable morph transition.
    if doTransitionPlot && wtSampIndex < prevWtSampIndex && n > sampsPerPeriod
        if mod(transitionPlotIndex, transitionPlotRateLimit) == 0
            figure(fig2),
            plot(y(n - floor(sampsPerPeriod) + 1: n), 'k.-'), ...
                ylim([-1.25, 1.25]), ...
%                 xlim([0, 100]), ...
                title('Wavetable transition'), ...
                drawnow limitrate;
        end
        transitionPlotIndex = transitionPlotIndex + 1;
    end
end

soundsc(y, Fs);
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
    subplot(313), ...
    spectrogram(y, 512, 64, 512, Fs, 'yaxis'), ...
    ylim([0, 22]);

if doTfPlot
    tfPlot(y, Fs, .01);
end

%%
% close all; 
% Fs = 44100;
% F0 = 1234;
% p = 100;
% t = linspace(0, p*(1/F0), p*(Fs/F0))';
% x = sin(2 * pi * F0 * t);
% plot(t, x)
% X = abs(fft(x, Fs));
% figure; plot(linspace(0, Fs-1, Fs), X);

%%
% clear; close all;
% Fs = 44100;
%  t = 0:1/Fs:0.1-1/Fs;
%  x = cos(2*pi*10*t);
%  xdft = fft(x);
%  xdft = xdft(1:length(x)/2+1);
%  df = Fs/length(x);
%  freqvec = 0:df:Fs/2;
%  stem(freqvec,abs(xdft),'markerfacecolor',[0 0 1])

%%
function y = interpolateZeroth(x, readIndex)
    n = floor(readIndex);
    readIndices = wrapIndices(n, length(x));
    y = x(readIndices(1));
end
 
 %%
 function y = interpolateLinear(x, readIndex)
    % Get the indices of the samples we need from the input signal.
    % Start by getting the index of the ith sample.
    n = floor(readIndex);
    % Create a vector of the two sample indices we need for linear interpolation.
    readIndices = wrapIndices([n, n+1], length(x));
    
    % Get the fractional part of the read index.
    alpha = mod(readIndex, 1);
    % Calculate the sample amplitude coefficients.
    % coeffs(1) => nth sample
    % coeffs(2) => (n+1)th sample
    coeffs = [1 - alpha, alpha];

    y = coeffs * x(readIndices);
 end
 
 %%
 function y = interpolateCubic(x, readIndex)
    % Get the indices of the samples we need from the input signal.
    % Start by getting the index of the ith sample.
    n = floor(readIndex);
    % Create a vector of the four sample indices we need for cubic 
    % interpolation. Wrap as required.
    readIndices = wrapIndices([n-1, n, n+1, n+2], length(x));

    % Get the fractional part of the read index.
    alpha = mod(readIndex, 1);
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

%     y3(n) = coeffs(1)*x(readIndices(1)) + ...
%         coeffs(2)*x(readIndices(2)) + ...
%         coeffs(3)*x(readIndices(3)) + ...
%         coeffs(4)*x(readIndices(4));
    % One-liner equivalent to the above:
    y = coeffs * x(readIndices);
 end
 
 %%
 function indices = wrapIndices(indices, bufferLength)
 %WRAPINDICES Wrap sample indices that fall outside of the length of a buffer.
    for i=1:length(indices)
        if indices(i) < 1
            indices(i) = bufferLength + indices(i);
        elseif indices(i) > bufferLength
            indices(i) = mod(indices(i), bufferLength);
        end
    end
 end

 %%
 function y = computeMipmap(x, Fs, F0)
 %COMPUTEMIPMAP Compute a wavetable mipmap bandlimited to Nyqvist.
    % E.g. for F0 = 44100/2048 = 21.5332, 8ve spacing (2 * F0), Fc = pi/2
    % E.g. for F1 = 2*F0, Fc = pi/4
    
    sampsPerPeriod = Fs/F0;
    Fc = ((sampsPerPeriod/length(x))/2)*.9;
    % Create windowed sinc LPF
    width = 2^11;
    h = hanning(width + 1) .* sinc(Fc*(-width/2:1:width/2)'); 
    % Normalize DC gain
    h = h./sum(h);
    H = fft(h, length(x));
    X = fft(x);
    
    y = ifft(X .* H);
    y = real(y);
 end
 
%%
function f=midi2hz(m)
    f= 440 * exp ((m-69) * log(2)/12);
end