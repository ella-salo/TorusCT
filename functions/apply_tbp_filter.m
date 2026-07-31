function rec_filtered = apply_tbp_filter(rec_tbp, alpha, s, M, target_res, mode)
% APPLY_TBP_FILTER Apply Tikhonov filter to TBP reconstruction.
%
% Computes the filtered torus backprojection f^{alpha,s} = f_TBP * p_{alpha,M}^s
% where p_{alpha,M}^s is the truncated Tikhonov filter with box size M.
% Supported modes:
%   'convolution' : computes the convolution as a numerical integral
%   'fourier'     : computes the convolution via FFT and iFFT

%%
switch mode
    % cFTBT: compute the numerical integral of the convolution
    case 'convolution'
        N = size(rec_tbp,1);
        dx = 1 / N;

        x_vec = linspace(0, 1-dx, N);
        [X1, X2] = meshgrid(x_vec, x_vec);

        % Build kernel p_alpha on grid
        p_kernel = zeros(N, N);

        for k1 = -M:M
            for k2 = -M:M
                k_norm_sq = 1 + k1^2 + k2^2;
                weight = 1 / (1 + alpha * k_norm_sq^s);

                p_kernel = p_kernel + weight * exp(2*pi*1i*(k1*X1 + k2*X2));
            end
        end

        p_kernel = real(p_kernel);

        % Convolution integral numerically
        rec_filtered = zeros(N, N);

        parfor i = 1:N
            for j = 1:N

                x1 = x_vec(i);
                x2 = x_vec(j);

                val = 0;

                for m = 1:N
                    for n = 1:N

                        y1 = x_vec(m);
                        y2 = x_vec(n);

                        % periodic shift x - y (mod 1)
                        dx1 = x1 - y1;
                        dx2 = x2 - y2;

                        dx1 = dx1 - round(dx1);
                        dx2 = dx2 - round(dx2);

                        % kernel index (nearest grid point)
                        idx1 = mod(round(dx1*N), N) + 1;
                        idx2 = mod(round(dx2*N), N) + 1;

                        val = val + rec_tbp(m,n) * p_kernel(idx1, idx2);

                    end
                end

                rec_filtered(i,j) = val * dx^2;

            end
        end

case 'fourier'
    % fFTBP: apply filter in Fourier domain.
    % Compute FFT of f_TBP, multiply by filter coefficients, inverse FFT.

    F = fft2(rec_tbp);

    % Build filter matrix in FFT frequency ordering
    filter_mat = zeros(target_res, target_res);
    for k1_idx = 1:target_res
        for k2_idx = 1:target_res
            % Convert FFT index to signed frequency
            k1 = k1_idx - 1; 
            if k1 > target_res/2
               k1 = k1 - target_res; 
            end
            k2 = k2_idx - 1; 
            if k2 > target_res/2
                k2 = k2 - target_res;
            end

            % Apply filter only within box Z_M, zero outside
            if abs(k1) <= M && abs(k2) <= M
                k_norm_sq = 1 + k1^2 + k2^2;
                filter_mat(k1_idx, k2_idx) = 1 / (1 + alpha * k_norm_sq^s);
            end
        end
    end

    rec_filtered = real(ifft2(F .* filter_mat));


end