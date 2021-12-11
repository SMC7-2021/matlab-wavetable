function y = computeMipmap(x, Fs, oversample, F0, mipmapsPerOctave)
    %COMPUTEMIPMAP Compute a wavetable mipmap bandlimited to Nyqvist.
    % E.g. for F0 = 44100/2048 = 21.5332, 8ve spacing (2 * F0), Fc = pi/2
    % E.g. for F1 = 2*F0, Fc = pi/4

    % Number of samples in one period of a signal at F0.
    sampsPerPeriod = Fs/F0;
    % Ratio of mipmap basis frequencies.
    interMipmapFreqRatio = 2^((12/mipmapsPerOctave)/12);
    % Calculate cutoff frequency.
    % For fundamental wavetable, 2048 samples at 44100, F0 = 21.5332 Hz, Fc is 
    % cutoff for the highest intended frequency for the mipmap, i.e.
    %   Fmax = F0 * mipmapFreqRatio
    %   Fs = 2 * (pi rad/sample)
    %   For F0, Fc = 1 * (pi rad/sample), but since we need the cutoff for the highest
    %   intended frequency:
    %   Fc = 1 / mipmapFreqRatio * (pi rad/sample)

    % For higher basis frequencies, need to take into account undersampling,
    % i.e. the need to skip some samples in order to get the desired frequency.

    % Get ratio by which current mipmap's bottom frequency undercuts the true
    % wavetable length.
    mipmapToLtRatio = (sampsPerPeriod) / length(x);
    % Cutoff is 1 over the ratio of the top/bottom of the mipmap, multiplied by
    % the undercut factor. Multiply by an oversampling-adjusted scalar to leave 
    % some headroom.
    % If no oversampling, set cutoff very high (.95 * nyquist) -- this becomes
    % our anti-alias filter.
    headroom = .95;
    % If oversampling, we can afford to filter higher, as the decimator will
    % take care of anti-alias filtering.
    if oversample > 1
        headroom = 3/(2*oversample);
    end
    Fc = (mipmapToLtRatio / interMipmapFreqRatio) * headroom;
    % Create windowed sinc LPF: sin(Fc * x)/(Fc * x)
    width = length(x);
    h = hann(width + 1) .* sinc(Fc*(-width/2:1:width/2)');
    % Normalize the filter.
    h = h./sum(h);
    H = fft(h, length(x));
    X = fft(x);

    % Multiply the frequency domain signal by the frequency response of the
    % filter, and return to the time domain.
    y = ifft(X .* H);
    y = real(y);
    % For reasons relating to the cyclical nature of the DFT, y is a phase
    % shifted version of the filtered input wavetable. I don't know well enough
    % why this is, however, here's a quick-and-dirty fix.
    y = [y(length(y)/2 + 1:end); y(1:length(y)/2)];
end