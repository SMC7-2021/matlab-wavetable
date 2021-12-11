%% Alias-Suppressed Oscillators Based on Differentiated Polynomial Waveforms
% Välimäki et al.
clear, close all;
Fs = 1000;
duration = 1;
F0 = 5;
t = linspace(0, 1, Fs)';
x = zeros(Fs*duration, 10);

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
    end
    subplot(515), plot(t(5:end), x(5:end, 10)), title 'Differentiated';
    pause(1);
end

% Detail of the amplitude transition due to differentiation.
figure;
stem(t(175:225), x(175:225, 10)), title('Detail of differentiated signal');