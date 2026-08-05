function k_data_primitive = map_sinogram_to_torus_primitive(sino_unit, det_px_d, N, X, thetas, xx, yy)
% MAP_SINOGRAM_TO_TORUS_PRIMITIVE Map sinogram data to torus Fourier samples.
%
% For each primitive direction (v1,v2), constructs the corresponding
% geodesic on the torus and samples projection data along it.
% Axis-aligned directions are handled directly, while other directions
% use geometric distance calculations and interpolation. The result is
% a matrix with [k1, k2, Dv] where Dv is the torus projected data in the
% direction v perpendicular to k.
%%
% Signed distance of point to a line
d2line = @(k,c,xy) (k.*xy(:,1)-xy(:,2)+c)./sqrt(k^2+1);

k_data_primitive = [];

totalNum1 = length(xx);

% Loop through primitive geodesic directions v
for ii = 1:length(xx)

    if mod(ii,10) == 0
        disp(['Interpolating: ' num2str(ii), '/', num2str(totalNum1)])
    end

    v1 = xx(ii);
    v2 = yy(ii);

    k1 = -v2;
    k2 = v1;

    if v2 == 0
        % data along horizontal axis, interpolation distances are the
        % discretization points
        p = thetas == 0;
        if sum(p) ~= 1
            warning('sum(p) ~= 1 with v = (%i,%i)',v1,v2);
        end
        proj_vals = interpolateProjection(sino_unit(:,p),(X(:,1)-0.5),det_px_d);
        Xdata = proj_vals;
    elseif v1 == 0
        % data along vertical axis, interpolation distances are the
        % discretization points
        p = thetas == -90;
        if sum(p) ~= 1
            warning('sum(p) ~= 1 with v = (%i,%i)',v1,v2);
        end
        proj_vals = interpolateProjection(sino_unit(:,p),(X(:,1)-0.5),det_px_d);
        Xdata = proj_vals;
    else
        % other directions than axes, interpolation distances are
        % solved based on the constructed geodesic
        p = (xx(:) / v1) == (yy(:) / v2);
        if sum(p) ~= 1
            warning('sum(p) ~= 1 with v = (%i,%i)',v1,v2);
        end

        Xdata = zeros(size(X,1),1); % allocations
        for jj = 1:size(X,1) % for all starting points on x-axis
            v = [xx(p) yy(p)]; % direction
            s = X(jj,:); % starting point

            % copies of [0,1]^2 to fit the geodesic vector, and the
            % center coordinates of these squares
            [x_grid, y_grid] = xy_grid(v,s);
            px_centers = [x_grid(:) y_grid(:)] + 0.5;

            % solve the distances of centers to the line
            k = v(2)/v(1);
            c = s(2) - k*s(1);
            d = d2line(k,c,px_centers);
            pxc_d = [px_centers d];
            pxc_d = pxc_d(abs(d)<sqrt(1/2),:);

            % interpolate the projection data for these distances
            proj_vals = interpolateProjection(sino_unit(:,p),pxc_d(:,3),det_px_d);

            % save the torus data
            Xdata(jj) = sum(proj_vals)/(norm(v));

        end
    end
    k_data_primitive = [k_data_primitive; k1 k2 Xdata'];
end

% (0,0) term:
p = thetas == 0;
proj_vals = sino_unit(:,p);
s = abs(det_px_d(2)-det_px_d(1));
Xdata = ones(size(X,1),1) .* sum(proj_vals(:).*s);
k_data_primitive = [k_data_primitive; 0 0 Xdata'];

end
