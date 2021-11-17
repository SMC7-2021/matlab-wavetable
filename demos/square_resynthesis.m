%% Demonstate Fourier synthesis of one period of a square wave
clear, close all;

maxCoefficient = 2^4 - 1;
Fs = 10000;
length = 2^8;
F0 = Fs/length;
T0 = 1/F0;
omega = 2 * pi * F0;
t = linspace(0, T0, length)';

% x(t) = (0 +) 4/pi cos(wt - pi/2) + 4/3pi cos(3wt - pi/2) + ...
y = zeros(length, 1);

for k=1:2:maxCoefficient
    component = ((4 / (k * pi)) * cos(k * omega * t - pi/2));
    y = y + component;
    plot(t, y),
        ylim([-1.4, 1.4]),
        xlim([0, T0]),
        title(sprintf('F%d = %f rad/s', k, F0 * k));
    drawnow;
    pause(.01);
end

s = square(2 * pi * F0 * t);
plot(t, y),
    ylim([-1.4, 1.4]),
    xlim([0, ceil(T0)]),
    hold on,
%     grid on,
    plot(t, s),
%     plot(t, abs(y - s)),
    xlabel('radians'),
    ylabel('amplitude'),
    legend('Approximation to square wave', 'Ideal square wave', 'Error magnitude');
    
%%
maxCoefficient = 3;
y3 = zeros(length, 1);
for k=1:2:maxCoefficient
    component = ((4 / (k * pi)) * cos(k * omega * t - pi/2));
    y3 = y3 + component;
end

maxCoefficient = 15;
y15 = zeros(length, 1);
for k=1:2:maxCoefficient
    component = ((4 / (k * pi)) * cos(k * omega * t - pi/2));
    y15 = y15 + component;
end

sp(1) = subplot(221); ...
    plot(t, y3, t, s), ...
    ylim([-1.4, 1.4]), ...
    xlim([0, T0]), ...
    xlabel('time (s)'), ...
    ylabel('x_3(t), s(t)'), ...
    title('Square wave resynthesis, N = 3');
sp(2) = subplot(222); ...
    plot(t, y15, t, s), ...
    ylim([-1.4, 1.4]), ...
    xlim([0, T0]), ...
    xlabel('time (s)'), ...
    ylabel('x_{15}(t), s(t)'), ...
    title('Square wave resynthesis, N = 15');
sp(3) = subplot(223); ...
    plot(t, abs(y3 - s)), ...
    ylim([0, 1]), ...
    xlim([0, T0]), ...
    xlabel('time (s)'), ...
    ylabel('|x(t) - s(t)|'), ...
    title('Error magnitude');
sp(4) = subplot(224); ...
    plot(t, abs(y15 - s)), ...
    ylim([0, 1]), ...
    xlim([0, T0]), ...
    xlabel('time (s)'), ...
    ylabel('|x(t) - s(t)|'), ...
    title('Error magnitude');

% get(sp(1), 'Position')
% get(sp(2), 'Position')
% get(sp(3), 'Position')
% get(sp(4), 'Position')

set(sp(1), 'Position', [0.100, 0.5838, 0.37, 0.3412]);
set(sp(2), 'Position', [0.550, 0.5838, 0.37, 0.3412]);
set(sp(3), 'Position', [0.100, 0.11, 0.37, 0.3412]);
set(sp(4), 'Position', [0.550, 0.11, 0.37, 0.3412]);