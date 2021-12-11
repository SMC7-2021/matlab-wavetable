function y = makeBLIT(Fs, f0, durationS)
%MAKEBLIT Create and return a bandlimited impulse train.
    T = durationS; %sec
    N = T*Fs; %num Samples
    n = [1:N]'; %sample Index
    t = n/Fs; %time

    % period in Seconds
    T1 = 1./f0;
    %sampling Period
    Ts = 1./Fs;
    % (1/f0)*fs
    P = T1/Ts;
    % Number of Harmonics below Nyqvist
    M = 2*floor(P/2.)+1;
    
    y = makeSinc(n, P, M);
    tfPlot(y, 44100, 4*P/Fs);
end

