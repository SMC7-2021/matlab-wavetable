function corners = detectCorners(signal, threshold)
%DETECTCORNERS Detect 'corners' in an arbitrary signal
% i.e. discontinuities in the signal or its derivatives

    if nargin == 1
        threshold = .5;
    end
    % Max order of differentiation.
    order = 2;
    derivatives = zeros(length(signal), order + 1);
    derivatives(:, 1) = signal;
    
    for n=0:length(signal) + 1
        n_ = mod(n, length(signal));
        for o=2:order + 1
            if n_ == 0
                derivatives(end, o) = derivatives(end, o-1) - derivatives(end - 1, o-1);
            elseif n_ == 1
                derivatives(n_, o) = derivatives(n_, o-1) - derivatives(end, o-1);
            else
                derivatives(n_, o) = derivatives(n_, o-1) - derivatives(n_-1, o-1);
            end
        end
    end

    % Normalize the second derivative
    M = max(abs(derivatives(:, order + 1)));
    d = (1/M) * derivatives(:, order + 1);
    
    % Find values of the normalized second derivative that exceed the threshold.
    corners = abs(d) > threshold;
    
    figure; plot(derivatives(:, 1));
end

