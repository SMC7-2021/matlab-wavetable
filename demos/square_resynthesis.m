%% Demonstate Fourier synthesis of one period of a square wave
clear, close all;

numCoefficients = 2^4 - 1;
Fs = 2 * pi;
length = 2^8;
F0 = Fs/length;
T0 = 1/F0;
omega = 2 * pi * F0;
t = linspace(0, T0, length)';

% x(t) = (0 +) 4/pi cos(wt - pi/2) + 4/3pi cos(3wt - pi/2) + ...
y = zeros(length, 1);

for k=1:2:numCoefficients
    component = ((4 / (k * pi)) * cos(k * omega * t - pi/2));
    y = y + component;
    plot(t, y),
        ylim([-1.4, 1.4]),
        title(sprintf('F%d = %f rad/s', k, F0 * k));
    drawnow;
    pause(.01);
end

s = square(2 * pi * F0 * t);
plot(t, y),
    ylim([-1.4, 1.4]),
    hold on,
%     grid on,
    plot(t, s),
    plot(t, abs(y - s)),
    xlabel('radians'),
    ylabel('amplitude'),
    legend('Approximation to square wave', 'Ideal square wave', 'Error magnitude');