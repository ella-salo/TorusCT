function fhat = compute_fourier_coefficients(k_sino_data, dN)
% COMPUTE_FOURIER_COEFFICIENTS Compute Fourier coefficients from sinogram data.
%
% For each frequency pair (k1,k2), evaluates the corresponding coefficient
% using a 1D DFT along the appropriate direction. Precomputed exponentials
% are used for efficiency.

%%
% values of exp(-2*pi*i*k*l/N) stored into memory
exp_l = meshgrid(0:dN-1);
exp_k = meshgrid(0:dN-1)';
% first coordinate k, second coordinate l
exp_kl = exp(2*pi*1i*exp_k.*exp_l/dN);
iexp_kl = exp_kl.^(-1);

% allocate fhat
fhat = zeros(size(k_sino_data,1),1);

for i = 1:size(k_sino_data,1)
    k1 = k_sino_data(i,1);
    k2 = k_sino_data(i,2);
    Dv = k_sino_data(i,3:end);
    if k1 == 0 && k2 == 0
        fhat(i) = DFT(Dv,0,exp_kl);
    elseif k1 ~= 0
        if k1 > 0
            fhat(i) = DFT(Dv,k1,exp_kl);
        else
            fhat(i) = DFT(Dv,-k1,iexp_kl);
        end
    elseif k2 ~= 0
        if k2 > 0
            fhat(i) = DFT(Dv,k2,exp_kl);
        else
            fhat(i) = DFT(Dv,-k2,iexp_kl);
        end
    else
        error('Unexpected case: k1 = %i, k2 = %i', k1, k2);
    end
end

end

