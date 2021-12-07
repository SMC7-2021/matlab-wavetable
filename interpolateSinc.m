function y = interpolateSinc(x, readIndex)
    % Get the fractional part of the read index.
    alpha = mod(readIndex, 1);
    % Compute a windowed sinc function of arbitrary width, centred on the read
    % index.
    width = 2^5;
    indices = (-width/2 + 1:width/2);
    s = hann(width)' .* sinc(indices - alpha);

    % Get the indices of the samples we need from the input signal.
    % Start by getting the index of the ith sample.
    n = floor(readIndex);
    % Create a vector of the sample indices we need for sinc interpolation.
    % Wrap as required.
    readIndices = wrapIndices(indices + n, length(x));

    y = s * x(readIndices);
end