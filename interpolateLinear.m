function y = interpolateLinear(x, readIndex)
    % Get the indices of the samples we need from the input signal.
    % Start by getting the index of the ith sample.
    n = floor(readIndex);
    % Create a vector of the two sample indices we need for linear interpolation.
    readIndices = wrapIndices([n, n+1], length(x));

    % Get the fractional part of the read index.
    alpha = mod(readIndex, 1);
    % Calculate the sample amplitude coefficients.
    % coeffs(1) => nth sample
    % coeffs(2) => (n+1)th sample
    coeffs = [1 - alpha, alpha];

    y = coeffs * x(readIndices);
end