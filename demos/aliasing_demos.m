%% 
clear, close all;
Fs = 44100;
% F0 = -43700; % or 44500, or 400, or 87800, or 88600, or (n*Fs) ± F0, n ∈ Z
F0 = 400;
n = -3;
x = sin(2 * pi * ((n * Fs) + F0) * linspace(0, Fs, Fs)');
soundsc(x, Fs);

%% Show a delightful animated plot that illustrates aliasing
clear, close all;
Fs = 10;
t = linspace(0, 1, Fs + 1)';
t2 = linspace(0, 1, Fs * 100)';
increment = .015;
for f=increment:increment:20
    x = sin(2 * pi * f * t);
    x2 = sin(2 * pi * f * t2);
    
    phase = 0;
    if mod(f, Fs/2) < increment
        fa = 0;
    elseif f > Fs / 2 && f <= Fs
        fa = Fs - f;
        phase = -pi;
    elseif f > Fs && f <= 3*Fs/2
        fa = -Fs + f;
    elseif f > 3*Fs/2 && f <= 2*Fs
        fa = 2*Fs - f;
        phase = -pi;
    else
        fa = f;
    end
%     xAlias = max(abs(x)) * sin(2 * pi * fa * t2 + phase);
    xAlias = sin(2 * pi * fa * t2 + phase);
    
    plot(t2, x2), ...
        ylim([-1., 1.]), ...
        title(sprintf('F0 = %f Hz', f));
    hold on;
    ln = plot(t2, xAlias); ...
        ln.LineWidth = 2;
    stem(t, x);
    hold off;
    
    drawnow;
    pause(.01);
end

%% Audition and show the spectrogram of some aliasing artifacts.
clear, close all;
% Sample rate.
Fs = 44100;
% Duration of output.
duration = 6;
% Length of output in samples.
numSamps = Fs * duration;
% Number of harmonics to generate.
numHarmonics = 2;
% Set up time vector.
t = linspace(0, duration, numSamps)';
% Initialize output.
y = zeros(numSamps, 1);

% Generate and sum the harmonics.
for h=1:numHarmonics
    % Set up frequency as a vector swept from 0 to some multiple of the sample
    % rate.
    f = linspace(0, Fs * h, numSamps)';
    % Compute the signal.
    x = sin(2 * pi * f .* t);
    % Sum it with the existing output.
    y = y + x;
end

% Plot the spectrogram; note the aliases and folded aliases.
spectrogram(y, 512, 64, 512, Fs, 'yaxis');
soundsc(y, Fs);

