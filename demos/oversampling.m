%% Generate signal
clear, close all;

Fs = 100;
duration = 1;
numSamps = Fs * duration;
t = linspace(0, duration, Fs * duration)';
f = linspace(-numSamps/2, numSamps/2 - 1, numSamps);
F0 = 4;
F0_2 = 7;
L = F0_2;
M = F0;

x = square(2 * pi * F0 * t);

figure('Position', [500, 10, 1500, 1000]);
subplot(721), plot(t, x), title('Square wave'), xlabel('Time (s)'), ylabel('Amplitude');

X = fftshift(fft(x));
subplot(722), plot(f, abs(X)), xlabel('Frequency (Hz)'), ylabel('Magnitude');

%% Upsample (interleaved with duplicates)
tz = linspace(0, duration, Fs * duration * L);
xz = zeros(size(x, 1) * L, size(x, 2));
fz = linspace(-(numSamps*L)/2, (numSamps*L)/2 - 1, numSamps * L);

for i=1:length(x)
    for j=1:L
        xz(j + (i-1) * L) = x(i) * L;
    end
end

subplot(723), plot(tz, xz), title('Upsampled, samples duplicated'), xlabel('Time (s)'), ylabel('Amplitude');

XZ = fftshift(fft(xz));
subplot(724), plot(fz, abs(XZ)), xlabel('Frequency (Hz)'), ylabel('Magnitude');

%% Upsample (interleaved with zeroes)
tz = linspace(0, duration, Fs * duration * L);
xz = zeros(size(x, 1) * L, size(x, 2));
% fz = linspace(-(numSamps)/2, (numSamps)/2 - 1, numSamps * L);
fz = linspace(-(numSamps*L)/2, (numSamps*L)/2 - 1, numSamps * L);

for i=1:length(x)
    xz(1 + (i-1) * L) = x(i) * L;
end
xz = upsample(x, L);

subplot(725), plot(tz, xz), title('Upsampled, interleaved with zeros'), xlabel('Time (s)'), ylabel('Amplitude');

XZ = fftshift(fft(xz));
subplot(726), plot(fz, abs(XZ)), xlabel('Frequency (Hz)'), ylabel('Magnitude');

%% Filter the upsampled version
xx = zeros(size(xz));

% Most basic FIR low pass filter.
for n=1:length(xz)
    if n > 1
        xx(n) = xz(n) + xz(n - 1);
    else
        xx(n) = xz(n);
    end
end

% for n=1:length(xz)
%     if n > 3
%         xx(n) = xz(n) + xz(n - 1) + xz(n - 2) + xz(n - 3);
%     elseif n > 2
%         xx(n) = xz(n) + xz(n - 1) + xz(n - 2);
%     elseif n > 1
%         xx(n) = xz(n) + xz(n - 1);
%     else
%         xx(n) = xz(n);
%     end
% end

% Slightly more sophisticated IIR filter.
% xx = filter([1,1], [1, -.5], xz);

subplot(727), plot(tz, xx), title('Basic first-order FIR low-pass'), xlabel('Time (s)'), ylabel('Amplitude');
XX = fftshift(fft(xx));
subplot(728), plot(fz, abs(XX)), xlabel('Frequency (Hz)'), ylabel('Magnitude');

%% Brick-wall frequency domain filter & downsample
XZ = fft(xz);
% figure;
% for i=length(x)/4:2*(length(x)*L)/4 - 1
for i=floor(length(x)*(M/L/2)):(length(x)*L/2) - length(x)*(M/L/2) - 1
    XZ(2*i) = 0.;
    XZ(2*i + 1) = 0.;
%     plot(abs(fftshift(XZ)));
%     drawnow;
end
xz_ = real(ifft(XZ));
y = downsample(xz_, M);

ty = linspace(0, duration * M/L, length(y))';
fy = linspace(-length(y)/2, length(y)/2 - 1, length(y))';
% figure('Position', [500, 10, 1500, 1000]);
subplot(729), plot(ty, y), title('Ablated in frequency domain'), xlabel('Time (s)'), ylabel('Amplitude');
Y = fftshift(fft(y));
subplot(7,2,10), plot(fy, abs(Y)), xlabel('Frequency (Hz)'), ylabel('Magnitude');

%% Downsample -- filter during decimation.
y2 = decimate(xz, M, 10);
% y2 = downsample(xz, M);
subplot(7,2,11), plot(ty, y2), title('Decimated with 10th order Chebyshev'), xlabel('Time (s)'), ylabel('Amplitude');
Y2 = fftshift(fft(y2));
subplot(7,2,12), plot(fy, abs(Y2)), xlabel('Frequency (Hz)'), ylabel('Magnitude');

%% Butterworth instead
cutoff = 6;
[b, a] = butter(6, cutoff/(Fs/2), 'low');
% zplane(b, a);
xz = filter(b, a, xz);
y2 = downsample(xz, M);
subplot(7,2,13), plot(ty, y2), title('6th order Butterworth, then downsampled'), xlabel('Time (s)'), ylabel('Amplitude');
Y2 = fftshift(fft(y2));
subplot(7,2,14), plot(fy, abs(Y2)), xlabel('Frequency (Hz)'), ylabel('Magnitude');