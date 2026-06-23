function [ dft_f ] = DFT( f_vec, k, exp_kl_mat )
% DFT Compute the k-th discrete Fourier coefficient using precomputed exponentials.

s = f_vec * exp_kl_mat( sub2ind(size(exp_kl_mat),(k+1)*ones(1,numel(f_vec)),1:numel(f_vec)) )';

dft_f = s/numel(f_vec);

end

