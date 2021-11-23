function [L, M] = getResamplingFactors(fsBase, fsTarget)
%GETRESAMPLINGFACTORS quickly get a rubbish approximation to a sampling rate ratio.
% Should use something like Diophantine equations, but that'll have to wait.
tau = fsTarget / fsBase;
L = 1;
p = -1;
% Let Matlab work out a decent approximation to tau, as long as the upsampling
% factor isn't 1.
while L <= 1
    p = p-1;
    [L, M] = rat(tau, 10^p);
end

% % M = round(tau * 10);
% M = fsBase;
% % L = 10;
% L = round(fsTarget);
% 
% mFrac = rem(M / 2, 1);
% lFrac = rem(L / 2, 1);
% 
% while(mFrac < .1 && lFrac < .1)
%     M = M/2;
%     L = L/2;
%     
%     mFrac = rem(M / 2, 1);
%     lFrac = rem(L / 2, 1);
% end
% end

