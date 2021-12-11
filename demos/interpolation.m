clear, close all;

Fs = 1000;
duration = 5;
F0 = 10;
order = 200;
t = linspace(0, duration, Fs * duration)';
x = square(2 * pi * F0 * t);

for p=1:order
    y = zeros(length(x), 1);
    for n=p+1:length(x)-p
        y(n) = x(n);
        for m=1:p
            y(n) = y(n) + (-1^(p-1))*(1/sqrt(p + m))*(x(n-m) + x(n+m));
        end
    end
    
    plot(t, y), ...
        title(sprintf('Order = %d', p));
    drawnow;
end

%%
t = linspace(-25, 25, 501)';
s = sinc(t) .* hann(length(t));% kaiser(length(t), 10);
wvtool(s);
% plot(t, s);

%%
clear, close all;

t = linspace(0, 100, 100);
x = square(2 * pi * 1 * t)';
ts = linspace(-5,100,600);
[Ts,T] = ndgrid(ts,t);
y = sinc(Ts - T)*x;

plot(t,x,'o',ts,y)
xlabel Time, ylabel Signal
legend('Sampled','Interpolated','Location','SouthWest')
legend boxoff

%%
clear, close all;
Fs = 44100;
Fc = 15000;
F0 = 5;
t = linspace(0, 1, Fs)';
x = square(2 * pi * F0 * t);
plot(t, x);
s = sinc(linspace(-2*pi, 2*pi, 100)') .* hann(100);
ts = 2 * (Fc / Fs) * ([1:50] - 50/2);
s1 = sinc(ts') .* hann(length(ts));
% 2 * Fc / Fs * (np.arange(N) - (N - 1) / 2)

% Create over-sampled sinc. Use that to interpolate wt.
wtLength = 2048;
sincLength = wtLength * 2^3;
nWt = linspace(0, (2 * pi) - 1/wtLength, wtLength)';
wt = square(nWt);
% x1 = audioread('./wavetables/violin_wt.wav');
% wt = resample(x1, wtLength, length(x1));
nSinc = linspace(-2*pi, 2*pi, sincLength)';
s = sinc(nSinc) .* hann(sincLength);
plot(nWt, wt, nSinc, s);

Fs = 44100;
F0 = 1000;
duration = 2;
sampsPerPeriod = Fs / F0;
wtStepsPerSample = wtLength / sampsPerPeriod;
y = zeros(Fs * duration, 1);
wtSampIndex = 1;
sincRange = 50;
sincWidth = .5;

tic

for n=1:Fs * duration
    for i=-sincRange/2:sincRange/2
        wtIndex = mod(round(wtSampIndex) + (i * round((wtLength*sincWidth) / sincRange)), length(wt)) + 1;
        sincIndex = floor((length(s) / 2) + (i/sincRange)*(length(s) - 1) + 1);
        y(n) = y(n) + (wt(wtIndex) * s(sincIndex));
    end
    
    wtSampIndex = mod(wtSampIndex + wtStepsPerSample, wtLength);
end

toc

figure;
plot(y(1:256));
soundsc(y, Fs);
figure;
plot(linspace(0, Fs, length(y))', abs(fft(y)));
figure;
spectrogram(y, 512, 64, 512, Fs, 'yaxis');


%% Alias-Suppressed Oscillators Based on Differentiated Polynomial Waveforms
% Välimäki et al.
clear, close all;
Fs = 1000;
duration = 1;
F0 = 5;
t = linspace(0, 1, Fs)';
x = zeros(Fs*duration, 10);

x(:, 1) = 2 * mod(linspace(0, duration * F0, Fs)', 1) - 1;
x(:, 1) = [repmat(audioread('./wavetables/vox_wt.wav'), 5, 1); zeros(Fs - 915, 1)];
% x(:, 1) = sin(2 * pi * F0 * t);
x(:, 1) = sawtooth(2 * pi * F0 * t, 1/2);
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

% and again
% x(:, 5) = (x(:, 1).^5)/120 - (x(:, 3).^3)/36 + x(:, 1)*(7/360);
% subplot(515), plot(t, x(:, 5)), title Fifth-order;

% Differentiate back to the sawtooth
x(:, 10) = x(:, 4);
for order=0:2
    x(:, 9) = x(:, 10);
    x(1, 10) = 0;
    for n=2:length(x)
        x(n, 10) = x(n, 9) - x(n-1, 9);
    end
    subplot(515), plot(t(5:end), x(5:end, 10)), title 'Differentiated';
    pause(.5);
end

% % Experiment with FIR filtering the sawtooth to smooth it out.
% x(:, 9) = x(:, 1);
% order = 2;
% for i=1:20
%     z = x(:, 9);
%     for n=order + 1:length(x)
%         x(n, 9) = (z(n) + z(n-1)) / 2;
%     end
% end
% figure, plot(x(:, 9));

%%