function y = runningAvg(x, alpha)
%RUNNINGAVG IIR aproximation of a running average
    if nargin == 1
        alpha = .995;
    end
    y = filter([1-alpha,0],[1,-alpha],x);
end

