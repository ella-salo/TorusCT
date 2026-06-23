function [x_grid, y_grid] = xy_grid(v,s)
% XY_GRID Construct grid of unit squares intersected by a geodesic.
%
% Given direction v and starting point s (on an axis), generates the grid
% of pixel indices covering one period of the discrete geodesic. The grid
% depends on the slope and orientation of the direction vector.
%%
if s(1) == 0
    % starting point at y-axis
    if sign(v(1)) > 0 && sign(v(2)) > 0
        % increasing slope
        if s(2) >= 0
            [x_grid, y_grid] = meshgrid( 0:v(1)-1, 0:v(2)+ceil(s(2))-1);
        else % s(1) < 0
            [x_grid, y_grid] = meshgrid( 0:v(1)-1, -1:v(2)-1);
        end
    elseif sign(v(1)) > 0 && sign(v(2)) < 0
        % decreasing slope
        if s(2) >= 0
            [x_grid, y_grid] = meshgrid( 0:v(1)-1, v(2)-1:0);
        else % s(1) < 0
            [x_grid, y_grid] = meshgrid( 0:v(1)-1, v(2)-1:0);
        end
    else
        error('slope not implemented');
    end 
elseif s(2) == 0
    % starting point at x-axis   
    if sign(v(1)) > 0 && sign(v(2)) > 0
        % increasing slope
        if s(1) >= 0
            [x_grid, y_grid] = meshgrid( 0:v(1)+ceil(s(1)) -1, 0:v(2)-1);
        else % s(1) < 0
            [x_grid, y_grid] = meshgrid( floor(s(1)):v(1)+floor(s(1)) , 0:v(2)-1);
        end
    elseif sign(v(1)) > 0 && sign(v(2)) < 0
        % decreasing slope
        if s(1) >= 0
            [x_grid, y_grid] = meshgrid( 0:v(1)+ceil(s(1)) -1, v(2):-1);
        else % s(1) < 0
            [x_grid, y_grid] = meshgrid( floor(s(1)):floor(s(1))+v(1), v(2):-1);
        end
    else
        error('slope not implemented');
    end    
else
    error('Only starting points along axes implemented.');
end

end
