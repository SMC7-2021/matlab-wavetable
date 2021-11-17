function [L, M] = getResamplingFactors(fsBase, fsTarget)
%GETRESAMPLINGFACTORS quickly get a rubbish approximation to a sampling rate ratio.
% Should use something like Diophantine equations, but that'll have to wait.
tau = fsTarget / fsBase;

% M = round(tau * 10);
M = fsBase;
% L = 10;
L = round(fsTarget);

mFrac = rem(M / 2, 1);
lFrac = rem(L / 2, 1);

while(mFrac == 0 && lFrac == 0)
    M = M/2;
    L = L/2;
    
    mFrac = rem(M / 2, 1);
    lFrac = rem(L / 2, 1);
end
end

