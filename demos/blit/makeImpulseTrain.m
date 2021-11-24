function y = makeImpulseTrain(Fs, f0, durationSecs)
%IMPULSETRAIN Generate an impulse train.
    % Number of samples to generate.
    N = durationSecs * Fs;

    % Time indices.
%     t = n/Fs;

    % Sample period.
    P = 1/f0 * Fs;
    
    % Output vector.
    y = zeros(ceil(N), 1);
    
    for n=1:N
        if floor(mod(n-1, P)) == 0
            y(n) = 1;
        else
            y(n) = 0;
        end
    end
    
    % Plot the (start of the) impulse train.
    tfPlot(y, Fs, 4*P/Fs);
end

