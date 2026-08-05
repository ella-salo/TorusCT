
%% Torus CT reconstruction from simulated data
% 26.6.2026

close all; clearvars; clc;

%% Configuration 
addpath('functions/')

% Method
% 'TorusCT'    : Fourier series in box of size (2N)^2
% 'StarTCT'    : Fourier series in extended star-shaped set K_{N,N2}
% 'TBP'        : Torus backprojection
method = 'StarTCT';

% Data
sample_mode = 'flag'; % Either 'shepp-logan', 'flag' or 'flag_rot30'
target_res = 256;
sample_res = 512;
noiselevel = 0.02; % in [0, 1]

% Fourier coefficient box size
N = 25;  % Fourier coefficients inside a box of size (2N)^2
N2 = 50; % StarTCT: extended frequency set K_{N,N2}

% Regularization parameters
alpha = 5*10^-6;
s = 1.25;

% TBP regularization
TBP_reg_method = 'fourier'; % Either 'convolution' or 'fourier'
M_filter = 100; % FTBP convolution filter size

% Positivity constraint
use_positivity = false;
% use_positivity = true;


%% Construct closed geodesic angles

[thetas, xx, yy] = create_closed_geodesic_angles(N);

% thetas \in [-90, 90)

figure
subplot(1,2,1)
hold on
[xx_all,yy_all] = meshgrid(-N:N,-N:N);
scatter(xx_all,yy_all,'k.');
scatter(xx,yy,'ko')
axis square
grid on
sgtitle('All coefficients and primitive rational directions')
subplot(1,2,2)
hold on
for i = 1:length(thetas)
        plot([0, cosd(thetas(i))], [0, sind(thetas(i))], 'k-')
end
axis equal
axis off

%% Generate the simulated data

[target, sample] = create_sample(sample_mode, target_res, sample_res);

[sino_for_torusCT, det_px_d] = radon(sample,thetas+90); % radon parametrizes angles from y-axis
sino_for_torusCT = sino_for_torusCT + noiselevel * max(abs(sino_for_torusCT(:))) .* randn(size(sino_for_torusCT));

%% RECONSTRUCTION ON THE TORUS

% Divide the sinogram by the resolution so that it represents measurements in [0,1]^2
sino_unit = sino_for_torusCT./sample_res;

% Detector pixel distances scaled to the unit square
det_px_d = (det_px_d/sample_res);

% Starting points of the geodesics
dN = 256;
X = [0:1/dN:1-1/dN ; zeros(1,dN)]';
Y = [ zeros(1,dN) ; 0:1/dN:1-1/dN]';

%% Map sinogram to torus data
% Solve the geodesics on the torus and interpolate corresponding data
% from the projections. Interpolate in the sense that the data is weighted
% sum of two closest projection pixels.

k_data_primitive = map_sinogram_to_torus_primitive(sino_unit, det_px_d, N, X, thetas, xx, yy);
mean_row = k_data_primitive(end, :);

%% Compute the Fourier coefficients and reconstruct

switch method

    case 'TorusCT'
        k_sino_data = [expand_k_data(k_data_primitive, 'box', N); mean_row];

        % Compute the Fourier coefficients
        fhat = compute_fourier_coefficients(k_sino_data, dN);

        % Build the Fourier series on the target_res grid
        [Xp,Yp] =  meshgrid(linspace(0,1,target_res),linspace(1,0,target_res));
        recf = FourierSeries(Xp,Yp,fhat,k_sino_data(:,1:2),alpha,s);

        % Correct for missing complex conjugates:
        % don't duplicate k=(0,0) term
        mean_ind = find(k_sino_data(:,1)==0 & k_sino_data(:,2)==0);
        rec = 2*real(recf)-k_sino_data(mean_ind,3);

    case 'StarTCT'
        
        k_sino_data = [expand_k_data(k_data_primitive, 'box', N2); mean_row];

        % Compute the Fourier coefficients
        fhat = compute_fourier_coefficients(k_sino_data, dN);

        % Build the Fourier series on the target_res grid
        [Xp,Yp] =  meshgrid(linspace(0,1,target_res),linspace(1,0,target_res));
        recf = FourierSeries(Xp,Yp,fhat,k_sino_data(:,1:2),alpha,s);

        % Correct for missing complex conjugates:
        % don't duplicate k=(0,0) term
        mean_ind = find(k_sino_data(:,1)==0 & k_sino_data(:,2)==0);
        rec = 2*real(recf)-k_sino_data(mean_ind,3);

    case 'TBP'
        % Torus backprojection (Theorem 2.3)
        % Each pixel value is obtained by identifying which geodesic in
        % direction v = k^perp passes through the reconstruction point and
        % reading off the precomputed torus data by linear interpolation.
        [Xp, Yp] = meshgrid(linspace(0,1,target_res),linspace(1,0,target_res));
        rec = zeros(target_res, target_res);
        n_terms = size(k_data_primitive, 1) - 1; % excluding (0,0) row

        for idx = 1:size(k_data_primitive, 1)
            k1_s   = k_data_primitive(idx, 1);
            k2_s   = k_data_primitive(idx, 2);
            g_vals = k_data_primitive(idx, 3:end);

            x_starts_ext = [X(:,1)'-1, X(:,1)', X(:,1)'+1];
            g_vals_ext   = [g_vals, g_vals, g_vals];

            % Geodesic direction is v = k^perp = (k2, -k1)
            % Find the x-axis starting point of the geodesic in direction v
            % passing through reconstruction point (x1, x2)
            v1 = k2_s;
            v2 = -k1_s;
            if abs(v2) > 0
                s_j = mod(Xp - Yp*(v1/v2), 1);
            else
                s_j = mod(Yp, 1);
            end
            val = reshape(interp1(x_starts_ext, g_vals_ext, s_j(:), 'linear'), ...
                          target_res, target_res);

            % Per-term mean correction: subtract (n_terms-1)/n_terms * term mean
            % so that the sum retains exactly one mean value and noise does not accumulate
            rec = rec + val - ((n_terms-1)/n_terms) * mean(g_vals);
        end
        % Apply Tikhonov filter if alpha > 0
        if alpha > 0
            rec = apply_tbp_filter(rec, alpha, s, M_filter, target_res, TBP_reg_method);
        end
end

% Positivity
if use_positivity
    rec = max(rec, 0);
end
%% Errors and visualization

error1 = 100*sum(abs(rec(:)-target(:))) / sum(abs(target(:)));
error2 = 100*sqrt( sum( (rec(:) - target(:)).^2 ) ) / sqrt(sum(target(:).^2));
errorinf = 100*max( abs( rec(:)-target(:) ) ) /max(abs(target(:)));

figure
imagesc(target)
colorbar('southoutside')
colormap gray
axis square
axis off
title('GT')

figure
subplot(1,2,1)
imagesc(rec)
title([method, ', error2 = ', num2str(error2)])
colorbar('southoutside')
colormap gray
axis square
axis off

subplot(1,2,2)
imagesc(target - rec)
title(['Difference between GT and ', method])
colorbar('southoutside')
axis square
axis off


%% Comparison to FBP

%% FBP reconstruction using the torus optimized angles

% For FBP we make a sinogram by interpolating from a higher
% resolution to avoid inverse crime
sino_for_FBP = create_radon_data_no_crime(target, sample, noiselevel, thetas+90);

recnfbp = iradon(sino_for_FBP, thetas+90);
recnfbp = recnfbp(2:end-1,2:end-1);

% Positivity
if use_positivity
    recnfbp = max(recnfbp, 0);
end

error1_fbp_torus_angles = 100*sum(abs(recnfbp(:)-target(:))) / sum(abs(target(:)));
error2_fbp_torus_angles = 100*sqrt(sum((recnfbp(:)-target(:)).^2)) / sqrt(sum(target(:).^2));
errorinf_fbp_torus_angles = 100*max(abs(recnfbp(:)-target(:))) / max(abs(target(:)));

% Visualize
figure
subplot(1,2,1)
imagesc(recnfbp)
colormap gray
axis square
axis off
colorbar('southoutside')
title(['FBP from Q_N, error2 = ', num2str(error2_fbp_torus_angles)])

subplot(1,2,2)
imagesc(target - recnfbp)
title('Difference between GT and FBP')
colormap gray
axis square
axis off
colorbar('southoutside')

%% FBP reconstruction using evenly distributed angles

N_ang = length(thetas); % Use same number of angles
even_thetas = -90 + [0:(N_ang-1)]/N_ang*180;

sino_equal = create_radon_data_no_crime(target, sample, noiselevel, even_thetas);

recnfbp_equal = iradon(sino_equal, even_thetas);
recnfbp_equal = recnfbp_equal(2:end-1,2:end-1);

% Positivity
if use_positivity
    recnfbp_equal = max(recnfbp_equal, 0);
end

error1_fbp_traditional = 100*sum(abs(recnfbp_equal(:)-target(:))) / sum(abs(target(:)));
error2_fbp_traditional = 100*sqrt(sum((recnfbp_equal(:)-target(:)).^2)) / sqrt(sum(target(:).^2));
errorinf_fbp_traditional = 100*max(abs(recnfbp_equal(:)-target(:))) / max(abs(target(:)));

figure
subplot(1,2,1)
imagesc(recnfbp_equal)
colormap gray
axis square
axis off
colorbar('southoutside')
title(['FBP, error2 = ', num2str(error2_fbp_traditional)])

subplot(1,2,2)
imagesc(target - recnfbp_equal)
title('Difference between GT and FBP')
colormap gray
axis square
axis off
colorbar('southoutside')


%% Error Table
M = [error1, error1_fbp_torus_angles, error1_fbp_traditional;
    error2, error2_fbp_torus_angles, error2_fbp_traditional;
    errorinf, errorinf_fbp_torus_angles, errorinf_fbp_traditional];

ErrorTable = array2table(M, ...
    'VariableNames', {method, 'FBP_torus_angles', 'FBP_traditional'}, ...
    'RowNames', {'L1_error','L2_error','Linf_error'});

disp(ErrorTable)
