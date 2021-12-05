%% Precalculated wavetables
% Variable length mipmaps. 
% NB, doesn't antialias as well as multiple_wavetable_7; some weirdness
% associated with jumps between the variable length mipmaps.

clear; close all; clc;
addpath('../helpers');
% Load some audio samples to use as wavetables.
x1 = audioread('./wavetables/vox_wt.wav');
x2 = audioread('./wavetables/violin_wt.wav');
% Sample rate.
Fs = 44100;
% Output duration.
outDurationS = .5;
% Output amplitude.
outAmp = 1.;
% Fundamental wavetable length.
Lt = 2^11;
shuffleWavetables = true;
% Number of mipmaps per octave.
mipmapDensity = 4;
% Sample interpolation type: 'truncate', 'linear', 'cubic'
interpolationType = 'cubic';
% Frequency output type: 'sweep10k', 'sweep5k', 'sweepMIDI', 'fixed440', 'fixed1k', 'fixed2k'
outputType = 'sweep5k';
% Max output samples to plot.
maxOutPlot = 1000;
% (Set up a rate-limiter for the wavetable transition plot.)
transitionPlotIndex = 0;
transitionPlotRateLimit = 100;
doTransitionPlot = true;
doTfPlot = false;

% Create a vector of output frequencies.
switch outputType
    case 'fixed440'
        F0 = linspace(440, 440, Fs * outDurationS)';
        doTfPlot = true;
    case 'fixed1k'
        F0 = linspace(1000.1, 1000.1, Fs * outDurationS)';
        doTfPlot = true;
    case 'fixed2k'
        F0 = linspace(2000, 2000, Fs * outDurationS)';
        doTfPlot = true;
    case 'sweep10k'
        outDurationS = 20;
        F0 = linspace(20, 10000, Fs * outDurationS)';
	case 'sweep5k'
        outDurationS = 20;
        F0 = linspace(20, 5000, Fs * outDurationS)';
    case 'sweepMIDI'
        outDurationS = 20;
        F0 = linspace(midi2hz(1), midi2hz(127), Fs * outDurationS)'; 
    otherwise % tweakable
        F0 = linspace(1500, 1500, Fs * outDurationS)';
        doTfPlot = true;
end

% Create an array of source wavetables. (Uncomment a statement to enable a
% wavetables; uncomment multiple statements to enable morphing.)
sourceWavetables = {
    % A sampled wavetable.
%     resample(x1, Lt, length(x1));
    % A sine wavetable.
%     sin(linspace(0, (2*pi) - (2*pi)/Lt, Lt)');
    % A discontinuous sine wavetable.
%     sin(linspace(0, 2.2*pi, Lt)');
    % A square wavetable.
    square(linspace(0, (2*pi) - (2*pi)/Lt, Lt)');
    % Another sampled wavetable.
%     resample(x2, Lt, length(x2));
    % A sawtooth wavetable.
    sawtooth(linspace(0, (2*pi) - (2*pi)/Lt, Lt)');
    % A noise wavetable
%     (rand(Lt, 1) * 2) - 1;
};
if shuffleWavetables
    sourceWavetables = sourceWavetables(randperm(length(sourceWavetables)));
end

% Plot the source wavetables.
figure('Name', 'Wavetables', 'Position', [100, 500, 600, 900]);

for i=1:length(sourceWavetables)
    subplot(length(sourceWavetables), 1, i), ...
        plot(1:Lt, sourceWavetables{i}, 'r-'), ...
        title(sprintf('Wavetable %d', i)), ...
        xlim([1, Lt]), ...
        ylim([-1.1, 1.1]), ...
        ylabel('amp.'), ...
        xlabel('sample index');
end

% Create a placeholder for the wavetable mipmaps.
wavetables = cell(1, length(sourceWavetables));

% Create a vector of frequencies to serve as the basis for the wavetable mipmaps.
% Start from the 'fundamental' frequency of a wavetable of length 'Lt' at
% sampling rate Fs (44100/2048 = 21.5332…) and proceed thnough nine octaves
% according to mipmap density.
mipmapFreqRatio = 2^((12/mipmapDensity)/12);
basisF0s = (Fs/Lt) * mipmapFreqRatio.^(0:9*mipmapDensity);
% Create a vector of simulated sampling rates for mipmap generation.
% For 'fundamental' wavetable frequency, use Fs. For higher mipmaps, use 
% basisFs = 

% Set up wavetable mipmaps.
for i=1:length(sourceWavetables)
    wt = sourceWavetables{i};
    x = cell(length(basisF0s), 1);
    
    % For each fundamental frequency, compute a mipmap of the wavetable 
    % bandlimited to Nyqvist.
    for f = 1:length(basisF0s)
        x{f} = computeMipmap(wt, Fs, basisF0s(f), mipmapDensity);
    end
    
    wavetables{i} = x;
end

% Placeholder for the wavetable morph samples.
wt = zeros(2, 1);
% To work with variable sample-length mipmaps, this demo defines phase as a 
% normalized quantity, 0 < phi < 1, and, to define the input sample indices, 
% multiplies phi by the length of the current mipmap.
phase = 0;
% Create output placeholder.
y = zeros(Fs * outDurationS, 1);

if doTransitionPlot
    fig2 = figure('Name', 'Wavetable transition', 'Position', [1000, 500, 600, 500]);
end

disp('Generating output...')

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
    % Default to mipmap 1 (for F0 below Fs/Lt)
    mipmap = 1;
    for f=length(basisF0s):-1:1
        if F0(n) > basisF0s(f)
            mipmap = f;
            break
        end
    end

    % Calculate the transitional wavetable samples.
    % First, convert normalized phase to read index.
    readIndex = (phase * length(wavetables{wtIndex}{mipmap})) + 1;
    switch interpolationType
        case 'truncate'
            wt(1) = interpolateZeroth(wavetables{wtIndex}{mipmap}, readIndex);
            wt(2) = interpolateZeroth(wavetables{nextWtIndex}{mipmap}, readIndex);
        case 'linear'
            wt(1) = interpolateLinear(wavetables{wtIndex}{mipmap}, readIndex);
            wt(2) = interpolateLinear(wavetables{nextWtIndex}{mipmap}, readIndex);
        otherwise % cubic
            wt(1) = interpolateCubic(wavetables{wtIndex}{mipmap}, readIndex);
            wt(2) = interpolateCubic(wavetables{nextWtIndex}{mipmap}, readIndex);
    end

    % Calculate the output sample.
    y(n) = outAmp * (...
        (1 - morph) * wt(1) + ...
        morph * wt(2) ...
    );

    % Calculate the number of samples per period of the wavetable to produce the
    % current frequency.
    sampsPerPeriod = Fs / F0(n);
    % Calculate the normalized phase increment.
    phaseIncrement = 1 / sampsPerPeriod;
    
    % (Previous index just used for plotting.)
    prevPhase = phase;
    
    % Upadate the wavetable sample index for the next iteration.
    phase = mod(phase + phaseIncrement, 1);

    % Make a cool plot of the wavetable morph transition.
    if doTransitionPlot && phase < prevPhase && n > sampsPerPeriod
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

disp('...done!')

soundsc(y, Fs);
% Plot beginning/end of output against time.
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

if doTfPlot
    tfPlot(y, Fs, .01);
end

%%
function y = interpolateZeroth(x, readIndex)
    % Truncate the read index.
    n = floor(readIndex);
    % Wrap as required.
    readIndices = wrapIndices(n, length(x));
    % Get output.
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
 function y = computeMipmap(x, Fs, F0, mipmapsPerOctave)
 %COMPUTEMIPMAP Compute a variable length wavetable mipmap bandlimited to Nyqvist.
    % E.g. for F0 = 44100/2048 = 21.5332, 8ve spacing (2 * F0), Fc = pi/2
    % E.g. for F1 = 2*F0, Fc = pi/4

    Lt = length(x);
    
    % Number of samples in one period of a signal at F0.
    sampsPerPeriod = Fs/F0;
    % Downsample input based on required level of detail. For fundamental
    % wavetable (F0 = 21.5332), downsample by factor of 1 (i.e. don't 
    % downsample); for F0*2, downsample by factor of 2; for F0*2^2, downsample
    % by factor of 2^2.
    periodsPerWT = Lt/sampsPerPeriod;
    powerTwo = floor(log2(periodsPerWT));
    M = 2^powerTwo;
    x = downsample(x, M);
    
    % Ratio of mipmap basis frequencies.
    mipmapFreqRatio = 2^((12/mipmapsPerOctave)/12);
    
    % Calculate cutoff frequency.
    % For fundamental wavetable, 2048 samples at 44100, F0 = 21.5332 Hz, Fc is 
    % cutoff for the highest intended frequency for the mipmap, i.e.
    %   Fmax = F0 * mipmapFreqRatio
    %   Fs = 2 * (pi rad/sample)
    %   For F0, Fc = 1 * (pi rad/sample), but since we need the cutoff for the highest
    %   intended frequency:
    %   Fc = 1 / mipmapFreqRatio * (pi rad/sample)
    
    % For higher basis frequencies, need to take into account undersampling,
    % i.e. the need to skip some samples in order to get the desired frequency.
    
    % Get amount by which current mipmap's bottom frequency is undersampled.
    % E.g. Lt = 2048, F0 = 25.6 Hz, therefore wavetable phase increment is 
    % greater than 1.
    underSamplingFactor = (sampsPerPeriod) / length(x);
    % Cutoff is 1 over the ratio of the top/bottom of the mipmap, multiplied by
    % the undersampling factor.
    Fc = (underSamplingFactor / mipmapFreqRatio);
    % Create windowed sinc LPF: sin(Fc * x)/(Fc * x)
    width = length(x);
    h = hann(width + 1) .* sinc(Fc*(-width/2:1:width/2)');
    % Normalize the filter.
    h = h./sum(h);
    H = fft(h, length(x));
    X = fft(x);
    
    % Multiply the frequency domain signal by the frequency response of the
    % filter, and return to the time domain.
    y = ifft(X .* H);
    y = real(y);
    % For reasons relating to the cyclical nature of the DFT, y is a phase
    % shifted version of the filtered input wavetable. I don't know well enough
    % why this is, however, here's a quick-and-dirty fix.
    y = [y(length(y)/2 + 1:end); y(1:length(y)/2)];
 end
