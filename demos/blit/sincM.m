function y = sincM(x, M)
%SINCM Bandlimited sinc function.
% y = sincM(x, M) Generate sinc of vector x with output limited to M harmonics.
    y = sin(pi*x) ./ (M*sin((pi*x)/M));
end

