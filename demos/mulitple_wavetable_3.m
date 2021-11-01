%% Multiple wavetable synthesis
% Now with more wavetables.
% Only supports fixed F0.

clear; close all;

% Constants
% Sample rate.
Fs = 44100;
% Output frequency.
F0 = 194.6;
% Output duration.
outDurationS = 9;
% Output amplitude.
outAmp = 1.;
% Wavetable length.
wtLength = 2^9;
% Max output samples to plot.
maxOutPlot = 1000;

% Load some audio samples.
x1 = audioread('./wavetables/vox_wt.wav');
x2 = audioread('./wavetables/violin_wt.wav');
x3 = audioread('./wavetables/drum_wt.wav');

% Create an array of wavetables.
wavetables = {
    % A sampled wavetable.
    resample(x1, wtLength, length(x1));
    % A sine wavetable
    sin(linspace(0, 2 * pi, wtLength)');
    % Another sampled wavetable.
    resample(x2, wtLength, length(x2));
    % A square wavetable.
    square(linspace(0, (2 * pi) - 1/wtLength, wtLength)');
    % Another sampled wavetable.
    resample(x3, wtLength, length(x3));
    % A sawtooth wavetable.
    sawtooth(linspace(0, (2 * pi) - 1/wtLength, wtLength)')
};
% Shuffle the wavetables for a laugh.
wavetables = wavetables(randperm(length(wavetables)));

% Plot them all.
figure('Name', 'Wavetables', 'Position', [100, 500, 600, 900]);

for i=1:length(wavetables)
    subplot(length(wavetables), 1, i), ...
        plot(1:wtLength, wavetables{i}, 'r-'), ...
        title(sprintf('Wavetable %d', i)), ...
        xlim([1, wtLength]), ...
        ylim([-1.1, 1.1]), ...
        ylabel('amp.'), ...
        xlabel('sample index');
end

fig2 = figure('Name', 'Wavetable transition', 'Position', [1000, 500, 600, 500]);

% Create output placeholder.
y = zeros(Fs * outDurationS, 1);

% Calculate the number of samples per period of the wavetable to produce the
% desired frequency.
sampsPerPeriod = Fs / F0;
wtStepsPerSample = wtLength / sampsPerPeriod;

% Placeholder for the transitional wavetable samples.
wt = zeros(2, 1);
% Initialize the wavetable sample index.
wtSampIndex = 1;
% (Set up a rate-limiter for the wavetable transition plot.)
transitionPlotIndex = 0;
transitionPlotRateLimit = 8;

% Transition linearly between the wavetables while writing to the output vector.
for n=1:(Fs * outDurationS)
    % Calculate the magnitudes for samples either side of the fractional index.
    magnitude2 = rem(wtSampIndex, 1);
    magnitude1 = 1 - magnitude2;
    
    % Wrap the wavetable sample index as necessary.
    prevSampIndex = floor(wtSampIndex);
    nextSampIndex = ceil(wtSampIndex);
    if prevSampIndex == 0
        prevSampIndex = wtLength;
    end
    if nextSampIndex > wtLength
        nextSampIndex = 1;
    end
    
    % Calculate how far through the output we are.
    outputDelta = n / (Fs * outDurationS);
    % Multiply this by the number of regions between wavetables (e.g. 5 
    % wavetables, 4 regions) to work out how far we are between the current 
    % wavetable pair.
    wtDelta = outputDelta * (length(wavetables) - 1);
    % The integer part is the index of the previous wavetable (+1 because matlab
    % isn't zero-indexed).
    wtIndex = floor(wtDelta) + 1;
    % The fractional part is the amount of transition between the pair of
    % wavetables.
    transition = rem(wtDelta, 1);
    
    % Calculate the index of the next wavetable.
    nextWtIndex = wtIndex + 1;
    % At the very end of output, wtDelta == length(wavetables) so nextWtIndex
    % will exceed the length of the wavetable array; so check for this.
    if nextWtIndex > length(wavetables)
        nextWtIndex = length(wavetables) - 1;
    end
    
    % Calculate the transitional wavetable samples.
    wt(1) = (1 - transition) * wavetables{wtIndex}(prevSampIndex) + ...
        transition * wavetables{nextWtIndex}(prevSampIndex);
    wt(2) = (1 - transition) * wavetables{wtIndex}(nextSampIndex) + ...
        transition * wavetables{nextWtIndex}(nextSampIndex);
    
    % Compose the output sample from summed, weighted transitional samples.
    y(n) = outAmp * ( ...
        magnitude1 * wt(1) + ...
        magnitude2 * wt(2) ...
    );
    
    % (Previous index just used for plotting.)
    prevWtSampIndex = wtSampIndex;
    
    % Upadate the wavetable sample index for the next iteration.
    wtSampIndex = mod(wtSampIndex + wtStepsPerSample, wtLength);
    
    % Make a cool plot of the wavetable transition.
    if wtSampIndex < prevWtSampIndex
        if mod(transitionPlotIndex, transitionPlotRateLimit) == 0
            figure(fig2),
            plot(y(n - floor(sampsPerPeriod) + 1: n), 'k.-'), ...
                ylim([-1, 1]), ...
                title('Wavetable transition'), ...
                drawnow limitrate;
        end
        transitionPlotIndex = transitionPlotIndex + 1;
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
