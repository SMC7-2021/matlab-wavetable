%% Multiple wavetable synthesis
% Supports frequency envelopes.
% Anti-aliasing windowed sinc filter applied.
% Aliasing reduced vs multiple_wavetable_4, but some other periodic artefacts
% introduced, possibly due to zeroth-order interpolation when determining the
% wavetable read index.

clear; close all;

% Constants
% Sample rate.
Fs = 44100;
% Output duration.
outDurationS = 9;
% Output amplitude.
outAmp = 1.;
% Wavetable length.
wtLength = 2^11;
% Max output samples to plot.
maxOutPlot = 1000;

% Load an audio sample.
x1 = audioread('./wavetables/drum_wt.wav');
x2 = audioread('./wavetables/violin_wt.wav');

% Create an array of wavetables.
wavetables = {
    % A sampled wavetable.
%     resample(x1, wtLength, length(x1));
    % Another sampled wavetable.
%     resample(x2, wtLength, length(x2));
    % A sine wavetable.
    sin(linspace(0, 2 * pi, wtLength)');
    % A square wavetable.
    square(linspace(0, (2 * pi) - 1/wtLength, wtLength)')
};

% Plot them.
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

% Create a vector of output frequencies.
F0 = linspace(66, 2027, Fs * outDurationS)';
% F0 = linspace(150, 350, Fs * outDurationS)' + sin(2 * pi * 1.5 * linspace(0, outDurationS, outDurationS * Fs)') .* ...
%     linspace(0, 30, Fs * outDurationS)';
F0 = linspace(21.5332, 2001.5332, Fs * outDurationS)';

% Placeholder for the transitional wavetable samples.
wt = zeros(2, 1);
% Initialize the wavetable sample index.
wtSampIndex = 1;
% The length of the sinc filter.
sincLength = wtLength * 2^3;
% Vector over which the sinc is defined.
nSinc = linspace(-2*pi, 2*pi, sincLength)';
% The windowed sinc function.
s = sinc(nSinc) .* hann(sincLength);
% Sinc filter range: the number of sinc coefficients to use for each sample.
sincRange = 50;
% (Set up a rate-limiter for the wavetable transition plot.)
transitionPlotIndex = 0;
transitionPlotRateLimit = 200;

tic

% Transition linearly between the wavetables while writing to the output vector.
for n=1:(Fs * outDurationS)
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
    if nextWtIndex == 0
        nextWtIndex = 1;
    end
    
    % Calculate the portion of the wavetable over which to apply the sinc
    % filter, based on the current frequency.
    % NB final multiplication is arbitrary -- need to revisit this.
    sincWidth = (F0(n) / (Fs/2)) * 8;
    
    % Apply the sinc filter by summing samples, centred on the 'current' sample,
    % weighted by values from the sinc function defined above.
    for i=-(sincRange/2):(sincRange/2)
        % Get the current wavetable read-index. 
        wtReadIndex = mod(...
            round(wtSampIndex) + (i * round((wtLength*sincWidth) / sincRange)), ...
            wtLength ...
        ) + 1;

        % Get the index from which to read the amplitude coefficient from the
        % sinc lookup.
        sincIndex = floor((sincLength / 2) + (i/sincRange)*(sincLength - 1) + 1);
        
        % Calculate the inter-wavetable sample.
        transitionalSample = (1 - transition) * wavetables{wtIndex}(wtReadIndex) + ...
            transition * wavetables{nextWtIndex}(wtReadIndex);
        
        y(n) = y(n) + (transitionalSample * s(sincIndex));
    end
    
    % Scale the output. Why 10/sincRange? More investigation required.
    y(n) = (10/sincRange) * y(n);

    % Calculate the number of samples per period of the wavetable to produce the
    % current frequency.
    sampsPerPeriod = Fs / F0(n);
    wtStepsPerSample = wtLength / sampsPerPeriod;
    
    % (Previous index just used for plotting.)
    prevWtSampIndex = wtSampIndex;
    
    % Upadate the wavetable sample index for the next iteration.
    wtSampIndex = mod(wtSampIndex + wtStepsPerSample, wtLength);

    % Make a cool plot of the wavetable transition.
    if wtSampIndex < prevWtSampIndex && n > sampsPerPeriod
        if mod(transitionPlotIndex, transitionPlotRateLimit) == 0
            figure(fig2),
            plot(y(n - floor(sampsPerPeriod) + 1: n), 'k.-'), ...
                ylim([-1, 1]), ...
%                 xlim([0, 100]), ...
                title('Wavetable transition'), ...
                drawnow limitrate;
        end
        transitionPlotIndex = transitionPlotIndex + 1;
    end
end

toc

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

tfPlot(y, Fs, .01);