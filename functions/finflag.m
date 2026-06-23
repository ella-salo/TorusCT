function [ v ] = finflag( xy, theta )
% FINFLAG
%   Returns gray value in a "nordic flag"

% Flag dimensions
% height 11 units
% width 18 units
% width of the cross 3 units
% height of the space below and above cross 4 units
% width of the inner space next to cross 5 units
% width of the outer space next to the cross 10 units

% allocate
v = zeros(size(xy,1),1);


% matrix for rotation
R = [cosd(theta)    -sind(theta)
     sind(theta)    cosd(theta)];
xy = (R*(xy' - 0.5*ones(size(xy')) ))' + 0.5*ones(size(xy));


% compute the gray values for all given coordinates
for i = 1:size(xy,1)
    x = xy(i,1);
    y = xy(i,2);

    if x < 0.14 || x > 0.86 || y < 0.28 || y > 0.72 % y < 0.28
        v1 = 0;
    else
        if y > 0.44 && y < 0.56
            v1 = 0.3;
        elseif x > 0.34 && x < 0.46
            v1 = 0.3;
        else
            v1 = 0.9;
        end
    end

    v(i) = v1;
    
end

% return
v = v';

end

