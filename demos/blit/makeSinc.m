function y = makeSinc(n, P, M)
%MAKESINC make a sampled sinc function

% y(n) = (M/P) * Sinc_M[(M/P)n]
% where P is the period in samples, and:
% Sinc_M(x) = sin(pi * x) / (M * sin(pi * x / M))
% M, number of harmonics:
% M = 2 * floor(P/2)+1 

    y = zeros(size(n));
    
    if M==P
        % If M==P, i.e. x in sincM is an integer, numerator of sincM will be 
        % zero. Only compute the inverse condition, i.e. after P has been 
        % adjusted (see below).
        cond = zeros(size(n));
        condInv = ones(size(n));
    else
        cond = mod(n / P, 1) ~= 0;
        % Denominator of sincM will be zero if P is an integer factor of n;
        % adjust P slightly for matching values of n.
        condInv = mod(n / P, 1) == 0;
    end
    
    % Do affirmative condtions.
    y(cond) = (M/P) * sincM((M/P)*n(cond), M);
    % Adjust P slightly.
    P = P - 1e-3;
    % Do negative conditions
    y(condInv) = (M/P) * sincM((M/P)*n(condInv), M);
end

