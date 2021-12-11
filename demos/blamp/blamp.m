clear, close all, clc;

addpath('../../helpers');
Fs = 44100;
duration = .1;
t = linspace(0, duration, duration*Fs)';

x = sawtooth(2 * pi * 400 * t, .5);
% x = sawtooth(2 * pi * 400 * t);
% x = square(2 * pi * 345.87 * t);
% x = sin(2 * pi * 200 * t);

% wt = audioread('../wavetables/vox_wt.wav');
% x = zeros(Fs*duration, 1);
% readIndex = 1;
% x = repmat(wt, 200, 1);
% for n=1:44100
%     alpha = mod(readIndex, 1);
%     nextReadIndex = floor(readIndex);
%     if nextReadIndex == 0
%         nextReadIndex = 1;
%     end
%     prevReadIndex = nextReadIndex-1;
%     if prevReadIndex == 0
%         prevReadIndex = length(wt);
%     end
%     x(n) = ((1-alpha) * wt(nextReadIndex)) + (alpha * wt(prevReadIndex));
%     readIndex = mod(readIndex + 15.1, length(wt));
% end

duration = length(x)/Fs;

corners = detectCorners(x, .5);
c = double(corners);
c(c==0) = NaN;
t = linspace(0, duration, duration*Fs)';
plot(t, c.*x, 'x', t, x), ...
    title('Corners'), ...
    grid on;

tfPlot(x, Fs, .01);
soundsc(x, Fs);
pause(duration);

y = polyBLAMP(x);
tfPlot(y, Fs, .01);
soundsc(y, Fs);