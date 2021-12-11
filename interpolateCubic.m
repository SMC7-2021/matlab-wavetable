function y = interpolateCubic(x, readIndex)
    % Get the indices of the samples we need from the input signal.
    % Start by getting the index of the ith sample.
    n = floor(readIndex);
    % Create a vector of the four sample indices we need for cubic 
    % interpolation. Wrap as required.
    readIndices = wrapIndices([n-1, n, n+1, n+2], length(x));

    % Get the fractional part of the read index.
    alpha = mod(readIndex, 1);
    % Calculate the sample amplitude coefficients.
    % coeffs(1) is the coefficient for the (n-1)th sample
    % coeffs(2) => nth sample
    % coeffs(3) => (n+1)th sample
    % coeffs(4) => (n+2)th sample
    coeffs = [ ...
        -alpha * (alpha - 1) * (alpha - 2) / 6, ...
        (alpha - 1) * (alpha + 1) * (alpha - 2) / 2, ...
        -alpha * (alpha + 1) * (alpha - 2) / 2, ...
        alpha * (alpha + 1) * (alpha - 1)/6 ...
    ];

    %     y = coeffs(1)*x(readIndices(1)) + ...
    %         coeffs(2)*x(readIndices(2)) + ...
    %         coeffs(3)*x(readIndices(3)) + ...
    %         coeffs(4)*x(readIndices(4));
    % One-liner equivalent to the above:
    y = coeffs * x(readIndices);
end