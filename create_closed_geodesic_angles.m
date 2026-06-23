function [thetas, xx, yy] = create_closed_geodesic_angles(N)
% CREATE_CLOSED_GEODESIC_ANGLES Generate coprime lattice directions (k1,k2)
% and corresponding projection angles, removing duplicate directions.

%%
[xx,yy] = meshgrid(0:N,-N:N);
% It is enough to compute the data in points whose gcd ~= 1 
gcd_xy = gcd(xx(:),yy(:));
gcd_xy_one = gcd_xy == 1;
xx( ~gcd_xy_one ) = [];
yy( ~gcd_xy_one ) = [];

% Remove the duplicate of (0,-1) which is (0,1)
xx(2) = [];
yy(2) = [];

% Projection directions based on points defined above
thetas = atand( yy ./ xx ) ;

end

