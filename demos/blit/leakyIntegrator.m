function y = leakyIntegrator(x, alpha)
%LEAKYINTEGRATOR IIR lowpass, 1 / (1 - a*z^-1) with alpha very close to 1
    if nargin == 1
        alpha = 1-10e-6;
    end
    
    y = filter([1.0, 0.], [1., -alpha], x);
end

