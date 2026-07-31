%% Torus CT reconstruction from real walnut data
% 26.6.2026

close all; clearvars; clc;

%% Configuration
addpath('functions/')

% Method
% 'TorusCT'    : Fourier series in box of size (2N)^2
% 'StarTCT'    : Fourier series in extended star-shaped set K_{N,N2}
% 'TBP'        : Torus backprojection
method = 'StarTCT';

% Fourier coefficient box size
N = 10;  % Fourier coefficients inside a box of size (2N)^2
N2 = 50; % StarTCT: extended frequency set K_{N,N2}

% Regularization parameters
alpha = 0;
s = 0;
% TBP regularization
TBP_reg_method = 'fourier'; % Either 'convolution' or 'fourier'
M = 100; % FTBP convolution filter size

% Positivity constraint
use_positivity = false;
% use_positivity = true;


%% Import data

load walnut_data/FullSizeSinograms.mat sinogram1200

target_res = 2296;
sample_res = 2296;


%% Preprocess the sinogram

sinogram1200double = double(sinogram1200);  % 2296 × 1200
sino_fan           = zeros(size(sinogram1200));

% Normalize the sinogram columnwise and take negative logarithm
numAngles = size(sinogram1200, 2);
for iii = 1 : numAngles

    bkgArea = sinogram1200double(1:100, iii);
    bkg = mean(bkgArea);

    sino_fan(:, iii) = -log(sinogram1200double(:, iii) / bkg);
end

% Plot the fan beam sinogram
figure
imagesc(sino_fan)
colormap(gray)
title('Fan beam sinogram')
colorbar('southoutside')


%% Define the parameters of the fan beam sinogram

% Source rotation increment (degrees)
ang_step = 0.3;

% Detector length (mm)
W = 114.8;

% Focus-origin-distance (mm)
FOD = 110;

% Focus-detector-distance (mm)
FDD = 300;

% Number of detectors
num_det = size(sino_fan,1);

% Detector pixel size (mm)
% pixelSize = 0.050;
pixelSize = W / num_det;

% Geometric magnification from fan-beam geometry
M = FDD/FOD;

% Effective pixel size due to magnification (mm)
pixelSizeEff = pixelSize / M;

% Ad hoc correction to sinogram to compensate for slightly misaligned
% center of rotation
sino_fan = circshift(sino_fan, 5);


%% Convert fan-beam projections to parallel beam projections

[sino_par, Ploc, Pangles] = fan2para(sino_fan, FOD, ...
    'FanSensorSpacing', pixelSizeEff, ...
    'ParallelSensorSpacing', pixelSizeEff, ...
    'ParallelRotationIncrement',ang_step,...
    'FanRotationIncrement', ang_step, ...
    'FanSensorGeometry', 'line',...
    'parallelCoverage','halfcycle');

sino = sino_par;
sino_equal = sino_par;


%% Construct closed geodesic angles

[thetas, xx, yy] = create_closed_geodesic_angles(N);

% thetas \in [-90, 90)


%% Select the closest perpendicular directions from the sinogram

sino_modified = zeros(size(sino,1),length(thetas));

% Pangles is the sinogram angle with respect to y-axis
% thetas is the geodesic segment angle with respect to x-axis

thetas_in_sino = mod(Pangles-90,180) - 90; % Wrap angle to interval [-90,90]

for ii = 1:length(thetas)

    % Search closest perpendicular angle from the sinogram
    [~, ind_of_theta] = min(abs(thetas_in_sino-thetas(ii)));

    % Change signed distances to be in correct order for the rays
    if Pangles(ind_of_theta) < 90
        sino_modified(:,ii) = sino(:, ind_of_theta);
    else
        sino_modified(:,ii) = sino(end:-1:1, ind_of_theta);
    end

end

sino = sino_modified;


% Now each column sino(:, thetas(i)) contains data which is perpendicular to
% geodesic segment at angle thetas(i) which we want for the torus
% reconstruction


%% FBP reconstruction using equispaced angles
% Compare with traditional equispaced angled measurements using FBP

recnfbp_equal = iradon(sino_equal, -Pangles, target_res);
target = recnfbp_equal;

figure
imagesc(recnfbp_equal)
colormap gray
axis square
axis off
colorbar('southoutside')
title(['FBP all angles, target'])


%% RECONSTRUCTION ON THE TORUS

% detector pixel distances
W_eff = W / M;
det_px_d = -Ploc./W_eff; % Flip direction to represent the setting

% Divide the sinogram by the resolution so that it represents measurements in [0,1]^2
sino_unit = sino./sample_res;

% Starting points of the geodesics
dN = 512;
X = [0:1/dN:1-1/dN ; zeros(1,dN)]';
Y = [ zeros(1,dN) ; 0:1/dN:1-1/dN]';


%% Map sinogram to torus data
% Solve the geodesics on the torus and interpolate corresponding data
% from the projections. Interpolate in the sense that the data is weighted
% sum of two closest projection pixels.

k_data_primitive = map_sinogram_to_torus_primitive(sino_unit, det_px_d, N, X, thetas, xx, yy);
mean_row = k_data_primitive(end, :);


%% Reconstruction


switch method

    case 'TorusCT'
        k_sino_data = [expand_k_data(k_data_primitive, 'box', N); mean_row];

        % Compute the Fourier coefficients
        fhat = compute_fourier_coefficients(k_sino_data, dN);
        k_sino_data = k_sino_data_torus;
        % Build the Fourier series on the target_res grid
        [Xp,Yp] = meshgrid(0:1/(target_res-1):1);
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
        [Xp,Yp] = meshgrid(0:1/(target_res-1):1);
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
        [Xp_grid, Yp_grid] = meshgrid(linspace(0,1,target_res));
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
            if abs(v1) > 1e-10
                s_j = mod(Yp_grid - Xp_grid*(v2/v1), 1);
            else
                s_j = mod(Xp_grid, 1);
            end
            val = reshape(interp1(x_starts_ext, g_vals_ext, s_j(:), 'linear'), ...
                          target_res, target_res);

            % Per-term mean correction: subtract (n_terms-1)/n_terms * term mean
            % so that the sum retains exactly one mean value and noise does not accumulate
            rec = rec + val - ((n_terms-1)/n_terms) * mean(g_vals);
        end
        % Apply Tikhonov filter if alpha > 0
        if alpha > 0
            rec = apply_tbp_filter(rec, alpha, s, M, target_res, TBP_reg_method);
        end
end

% Positivity
if use_positivity
    rec = max(rec, 0);
end


%% Errors and visualization

error1 = 100*sum(abs(rec(:)-target(:))) / sum(target(:));
error2 = 100*sqrt( sum( (rec(:) - target(:)).^2 ) ) / sqrt(sum(target(:).^2));
errorinf = 100*max( abs( rec(:)-target(:) ) ) /max(abs(target(:)));


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
title('Difference')
colorbar('southoutside')
axis square
axis off


%% Comparison to FBP

%% FBP reconstruction using the torus optimized angles

recnfbp = iradon(sino, -thetas, target_res);

% Positivity
if use_positivity
    recnfbp = max(recnfbp, 0);
end

error1_fbp_torus_angles = 100*sum(abs(recnfbp(:)-target(:))) / sum(target(:));
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


%% Error Table
M = [error1, error1_fbp_torus_angles;
    error2, error2_fbp_torus_angles;
    errorinf, errorinf_fbp_torus_angles];

ErrorTable = array2table(M, ...
    'VariableNames', {method, 'FBP_torus_angles'}, ...
    'RowNames', {'L1_error','L2_error','Linf_error'});

disp(ErrorTable)