%% Alias-Suppressed Oscillators Based on Differentiated Polynomial Waveforms
% Välimäki et al.
clear, close all;
Fs = 44100; %init 1000
duration = 1;
F0 = 2500; %init 5
t = linspace(0, 1, Fs)';
x = zeros(Fs*duration, 10);

%period for F0 in samples
P = Fs / F0;

x(:, 1) = 2 * mod(linspace(0, duration * F0, Fs)', 1) - 1;
% x(:, 1) = [repmat(audioread('./wavetables/vox_wt.wav'), 5, 1); zeros(Fs - 915, 1)];
% x(:, 1) = sin(2 * pi * F0 * t);
% x(:, 1) = sawtooth(2 * pi * F0 * t, 1/2);
subplot(511), plot(t, x(:, 1)), title 'Trivial Saw';

% Integrate x
x(:, 2) = (x(:, 1).^2)/2;
subplot(512), plot(t, x(:, 2)), title Parabolic;

% intergrate x again
x(:, 3) = (x(:, 1).^3 - x(:, 1))/6;
subplot(513), plot(t, x(:, 3)), title Cubic;

% and again
x(:, 4) = (x(:, 1).^4)/24 - (x(:, 1).^2)/12;
subplot(514), plot(t, x(:, 4)), title Fourth-order;

% Differentiate back to the sawtooth
x(:, 10) = x(:, 4);

for order=0:2
    x(:, 9) = x(:, 10);
    x(1, 10) = 0;
    
    for n=2:length(x)
        x(n, 10) = x(n, 9) - x(n-1, 9);
        
        %apply basic scaling factor
        x(n, 8) = x(n, 10) * (P^3 / 192);
        
        %apply improved scaling factor
        x(n, 7) = x(n,10) * 3*pi / (24 * (2*sin(pi/P)^3));
        
    end
    subplot(515), plot(t(5:end), x(5:end, 10)), title 'Differentiated';
    pause(1);
end

%plot the differetiated sawtooth with amplitude scaling
%NB - skips the first 4 samples (they look weird, but sound normal...??)
figure;
subplot(211), plot(t(4:end), x(4:end, 8)), title 'scaled: basic'
subplot(212), plot(t(4:end), x(4:end, 7)), title 'scaled: improved'

%plot spectrograms
figure;
spectrogram(x(:, 10), 512, 64, 512, Fs, 'yaxis');
title("after integration (before scaling)")
figure;
spectrogram(x(:, 1), 512, 64, 512, Fs, 'yaxis');
title("before integration")
figure;
spectrogram(x(:, 8), 512, 64, 512, Fs, 'yaxis');
title("8th matrix - basic scaling")
figure;
spectrogram(x(:, 7), 512, 64, 512, Fs, 'yaxis');
title("7th matrix - improved scaling")

% Detail of the amplitude transition due to differentiation.
figure;
stem(t(175:225), x(175:225, 10)), title('Detail of differentiated signal');

soundsc(x(:, 7), Fs);