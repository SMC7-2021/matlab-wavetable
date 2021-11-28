%% Bandlimited impulse train experiment.
% Essentially a port of https://github.com/znibbles/01-advanced-synthesis-techniques
clear; close all; clc;
addpath('../../helpers')

Fs = 44100;
useCumsum = true;

%% Create a sawtooth wave from an ideal impulse train.
x = makeImpulseTrain(Fs, 2790, .5);

% Get DC component.
avg = mean(x);
% Integrate (after subtracting DC) to create a simple sawtooth wave.
if useCumsum
    % cumsum() essentially applies a first order IIR lowpass.
    simpleSaw = cumsum(x - avg);
else
    simpleSaw = zeros(length(x), 1);
    simpleSaw(1) = x(1);
    for n=2:length(x)
        simpleSaw(n) = (x(n)-avg) + .99*simpleSaw(n-1);
    end
end

% Remove DC component of result.
simpleSaw = simpleSaw - mean(simpleSaw);

soundsc(simpleSaw, Fs);
tfPlot(simpleSaw, Fs, .005);

%% A brickwall filter in frequency domain is a sinc in time domain.
brickwall_filter = ones(1024, 1);
brickwall_filter(128:end) = 0;

figure;
H = real(ifft(brickwall_filter));
plot(H(1:128)), ...
    grid on;

%% So set up a bandlimited impulse train based on a sinc function...
% If our impulse is bandlimited, then by definition any signal we derive from it
% -- by integrating, i.e. filtering -- is bandlimited too.

%% Create a sawtooth wave from a band-limited impulse train
y = makeBLIT(Fs, 801, 1.);
% As above, integrate it to create a sawtooth wave
avg = mean(y);
if useCumsum
    blitSaw = cumsum(y - avg);
else
    blitSaw = zeros(length(x), 1);
    blitSaw(1) = y(1);
    for n=2:length(y)
        blitSaw(n) = (y(n)-avg) + (1-10e-4)*blitSaw(n-1);
    end
end

blitSaw = blitSaw - mean(blitSaw);

tfPlot(blitSaw, Fs, 1);
soundsc(blitSaw, Fs);

%% Leaky integrator
% All-in-all, works pretty well for f0 based on complex ratios, not for round 
% integers, e.g. 1200 Hz. (See second parameter to makeBLIT.)

% Get DC component of BLIT
newAvg = runningAvg(y, .99);
mean(newAvg)
% Integrate after removing the DC component
saw = leakyIntegrator(y - newAvg);
% Remove DC component of result
saw = saw - runningAvg(saw);

soundsc(saw, Fs);
tfPlot(saw, 44100);