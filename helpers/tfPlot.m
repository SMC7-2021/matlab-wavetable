function tfPlot(x, Fs, varargin)
%TFPLOT Create and display a time/frequency plot
    fontName = 'Times';
    fontSize = 12;

    N = length(x);
    t = (1:N)/Fs;
    
    figure('Position', [500, 300, 1000, 400]);
    
    subplot(121),...
        plot(t, x),...
        grid on,...
        xlabel('Time (sec)'), ...
        ylabel('Amplitude'), ...
        ylim([-1.5, 1.5]), ...
        set(gca,'fontsize',fontSize,'fontname',fontName);
    
    % Limit the x-range for the time plot if a limit (in seconds) has been 
    % specified.
    if nargin == 3 && isnumeric(varargin{1})
        xlim([0, varargin{1}]);
    end
    
    w = chebwin(length(x), 200);
    X = abs(fft(w.*x, Fs));
    X = db(X/max(X));
    subplot(122),...
%         semilogx(linspace(0, Fs-1, Fs), X), ...
        plot(linspace(0, Fs-1, Fs), X, 'LineWidth', 1), ...
        axis([0 Fs/2 -100 5]), ...
        xlabel('Freq (kHz)'),...
        ylabel('Magnitude (dB)'), ...
        % title('Spectrum'),...
        set(gca,'XTick',[0 5e3 10e3 15e3 20e3]), ...
        set(gca,'XTickLabel',[{'0'}, '5', '10' '15' '20']), ...
        set(gca,'YTick',-200:20:0), ...
        set(gca,'fontsize',fontSize,'fontname',fontName), ...
        grid on;
end

