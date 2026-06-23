function data = interpolateProjection(projection, distances, detector_pixels)
% INTERPOLATEPROJECTION Interpolate projection values at given distances.
%
% For each query point, finds the two nearest detector pixels and computes
% a linear interpolation based on distance weights. Values outside the
% detector range are set to zero (compact support)
%%
values = zeros(size(distances));
for i = 1:numel(distances)
    % Two closest detector pixels, idx(1) and idx(2)
    if distances(i) < min(detector_pixels) || distances(i) > max(detector_pixels)
        values(i) = 0; % compact support
    else
        [d, idx] = sort(abs(detector_pixels - distances(i)));

        % Values at the two closest detector pixels
        val1 = projection(idx(1));
        val2 = projection(idx(2));

        % Interpolate total value as weighted average of neighbour values,
        % where each weight is the proportional distance to the detector 
        % neighbour pixel
        w1 = d(2) / (d(1) + d(2));
        w2 = d(1) / (d(1) + d(2));

        values(i) = (w1*val1 + w2*val2); 
    end
end

data = values;

end