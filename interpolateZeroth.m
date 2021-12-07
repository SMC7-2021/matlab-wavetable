function y = interpolateZeroth(x, readIndex)
    % Truncate the read index.
    n = floor(readIndex);
    % Wrap as required.
    readIndices = wrapIndices(n, length(x));
    % Get output.
    y = x(readIndices(1));
end