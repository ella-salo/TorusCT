function [sino_FBP] = create_radon_data_no_crime(target, sample, noiselevel, thetas)
% CREATE_RADON_DATA_NO_CRIME Create synthetic Radon data without inverse crime.
%
% This implementation follows the approach of Müller and Siltanen,
% specifically the supplementary MATLAB script:
%   "XRA_NoCrimeData_comp.m"
% from the book:
%   Müller & Siltanen, Linear and Nonlinear Inverse Problems with Practical Applications.

%%
target_res = size(target,1);
res = size(sample,1);

% Radon transforms
[m_high, s_high] = radon(sample, thetas);
[m_low, s_low] = radon(target, thetas);

% Scaling and origin shift similarly to Mueller and Siltanen

ratio = target_res / res; % scaling factor
s_high_scaled = s_high * ratio; % scale coordinates

% Origin in low-res grid
x  = 0.5 + [0:target_res-1];
[X,Y] = meshgrid(x,x);
orind = floor((size(target)+1)/2);
orx  = X(orind(1),orind(2));
ory  = Y(orind(1),orind(2));

% Origin in high-res grid (scaled)
x2 = (0.5 + [0:res-1]) * ratio;
[X2,Y2] = meshgrid(x2,x2);
orind2 = floor((size(sample)+1)/2);
orx2 = X2(orind2(1),orind2(2));
ory2 = Y2(orind2(1),orind2(2));

% Shift
odist = sqrt((orx-orx2)^2 + (ory-ory2)^2);

% Interpolate high-res -> low-res sinogram (no crime)
m_interp = zeros(size(m_low));

for i = 1:length(thetas)
    shi = s_high_scaled(:) + odist * cosd(thetas(i) + 45);
    m_interp(:,i) = interp1(shi, m_high(:,i), s_low(:), 'pchip', 0);
end

% Magnitude correction
m_interp = m_interp * ratio;

% Add noise
sino_FBP = m_interp + noiselevel * max(abs(m_interp(:))) .* randn(size(m_interp));

end
