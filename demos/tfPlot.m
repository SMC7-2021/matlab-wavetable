function tfPlot(x, Fs, varargin)
%TFPLOT Create and display a time/frequency plot
    N = length(x);
    t = (1:N)/Fs;
    
    figure('Position', [500, 300, 1000, 400]);
    
    subplot(121),...
        plot(t, x),...
        grid on,...
        xlabel('Time[sec]');
    
    % Limit the x-range for the time plot if a limit (in seconds) has been 
    % specified.
    if nargin == 3 && isnumeric(varargin{1})
        xlim([0, varargin{1}]);
    end
    
    X = abs(fft(x));
    subplot(122),...
        semilogx(linspace(0, Fs/2, ceil(N/2)), X(1:ceil(N/2))), ...
        xlabel('Freq[Hz]'),...
        title('Spectrum'),...
        grid on;
end

