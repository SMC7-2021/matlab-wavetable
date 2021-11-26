Time = 0:0.001:5;
T = Time(2)-Time(1);             % Sampling period
Fs = 1/T;                        % Sampling frequency
L = numel(Time);                 % Length of signal
Disp_Data = 5*(sin(10*2*pi*Time)+0.1*((rand(1,L)-0.5)*2));
M_data = mean(Disp_Data);
RMS_data = rms(Disp_Data);
figure
plot(Time, Disp_Data)
hold on
plot(xlim, [1 1]*M_data, 'LineWidth',1.5)
plot(xlim, [1 1]*RMS_data, 'LineWidth',1.5)
hold off
grid
xlabel('Time')
ylabel('Amplitude')
legend('Signal', 'Signal Mean', 'Signal RMS')