function [value] = FourierSeries(X,Y,fhat,k,alpha,s)
% FOURIERSERIES Evaluate truncated Fourier series with smoothing weights.

value = 0; %initialize
for i = 1:size(k,1)
    
    P = (1 + alpha*(1 + norm(k(i,:))^2)^(s))^(-1);
    value = value + ...
        P*fhat(i)*exp(2*pi*1i*(k(i,1)*X + k(i,2)*Y));

end
