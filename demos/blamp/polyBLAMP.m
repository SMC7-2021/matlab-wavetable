function y = polyBLAMP(x, threshold)
%POLYBLAMP Apply four-coefficient polynomial band-limited ramp approximation to signal.

% Iterate over samples
% Check for a discontinuity in the first derivative
%   i.e. look for spikes in the second derivative
% Locate the sample where the discontinuity occurs
% Apply polyBLAMP centred on that sample
%   Check whether it's a rising or falling edge -- sign(secondDerivative)
%   Get magnitude of the discontinuity from the slope of the signal at the edge
%   Estimate where the discontinuity occurs (between samples)
%   Centre the BLAMP funciton on that fractional sample position

    % Allocate output vector
    y = zeros(size(x));

    % Control variables used to implemnt algorithm inside a loop
    d0 = zeros(2, 1);
    d1 = zeros(2, 1);
    d2 = 0;

    corner = 1;
    corner_n1 = 1;
    
    if nargin == 1
        threshold = .5;
    end

    flag = 0;

     xn = 0;
    xn1 = 0;
    xn2 = 0;
    xn3 = 0;

     yn = 0;
    yn1 = 0;
    yn2 = 0;
    yn3 = 0;

    % Main loop
    for n=1:length(x)

        xn = x(n);
        d0(1) = d0(2);
        d0(2) = xn;
        
        % Detect a corner
        % First derivative is current sample minus previous sample.
        d1(1) = d1(2);
        d1(2) = d0(2) - d0(1);
        % Second derivative is current difference minus previous difference.
        d2 = d1(2) - d1(1);

        if abs(d2) > threshold
            corner = 1;
        else
            corner = 0;
        end

        yn = xn;

        % If corner then implement correction
        if (corner-corner_n1) ~= 0
            flag = 1; 
        elseif flag==1
            flag = 0;

            % Determine whether the corner represents a rising or falling edge.
            pol = sign(xn1);
            % Determine the amplitude of the point of corner onset.
            L = xn3;

            % Perform polynomial interpolation around corner boundaries
            p_a =  (-1/6)*xn3  +    0.5*xn2  -  0.5*xn1  +  (1/6)*xn;
            p_b =         xn3  -  (5/2)*xn2  +    2*xn1  -    0.5*xn;
            p_c = (-11/6)*xn3  +      3*xn2  -  1.5*xn1  +  (1/3)*xn;
            p_e =         xn3;

            % Perform Newton-Raphson to estimate corner onset point
            x_d = 1.5;
%             for m=1:100
%                 err = (p_a*x_d^3 + p_b*x_d^2 + p_c*x_d + p_e - L) / ...
%                     (3*p_a*x_d^2 + 2*p_b*x_d + p_c);
%                 if abs(err) > 1e-6
%                     x_d = x_d - err;
%                 else
%                     break
%                 end
%             end

            % Estimate slope at corner
            mu = abs(3*p_a*x_d^2 + 2*p_b*x_d + p_c);

            % Fractional delay required to center polyBLAMP
            d = x_d - 1;

            % Compute polyBLAMP coefficients
            h0 = -d^5/120 + d^4/24 - d^3/12 + d^2/12 - d/24 + 1/120;
            h2 = -d^5/40 + d^4/24 + d^3/12 + d^2/12 + d/24 + 1/120;
            h1 = d^5/40 - d^4/12 + d^2/3 - d/2 + 7/30;
            h3 = d^5/120;

            % Superimpose polyBLAMP at corner
            yn3 = yn3 - pol*h0*mu;
            yn2 = yn2 - pol*h1*mu;
            yn1 = yn1 - pol*h2*mu;
             yn = yn  - pol*h3*mu;

        end

        % Output sample
        y(n) = yn3;

        % Update ALL control variables
        corner_n1 = corner;

        xn3 = xn2; xn2 = xn1; xn1 = xn;

        yn3 = yn2; yn2 = yn1; yn1 = yn;

    end
end

