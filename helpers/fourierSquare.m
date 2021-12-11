function y = fourierSquare(Fs, F0, duration)
%FOURIERSQUARE generate a perfectly bandlimited square wave approximation.
    samps = Fs * duration;
    nyquist = Fs/2;
    maxCoeff = nyquist/F0;
    t = linspace(0, duration, samps)';
    y = zeros(samps, 1);
    omega = 2 * pi * F0;
    
    for k=1:2:floor(maxCoeff)
        component = ((4 / (k * pi)) * cos(k * omega * t - pi/2));
        y = y + component;
    end
end