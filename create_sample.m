function [target, sample] = create_sample(mode, target_res, sample_res)
% CREATE_SAMPLE Generate high- and low-resolution test images.
%
% Supported modes:
%   'shepp-logan' : standard Shepp-Logan phantom
%   'flag'        : flag created by finflag.m
%   'flag_rot30'  : rotated version of the flag pattern (30 degrees)
%
%% Create high resolution sample
if strcmpi(mode, 'shepp-logan')
    sample = phantom(sample_res);

elseif strcmpi(mode, 'flag')
    flag_high = zeros(sample_res);
    X = 0:1/(sample_res-1):1;
    Y = X;
    rot = 0;
    for i = 1:numel(X)
        flag_high(:,i) = finflag( [ ones(size(X,2),1)*X(i)' Y' ], rot );
    end
    sample = flag_high;
elseif strcmpi(mode, 'flag_rot30')
    flag_high_rot = zeros(sample_res);
    X = 0:1/(sample_res-1):1;
    Y = X;
    rot = 30;
    for i = 1:numel(X)
        flag_high_rot(:,i) = finflag( [ ones(size(X,2),1)*X(i)' Y' ], rot );
    end
    sample = flag_high_rot;
end

%% Create low-resolution sample
if strcmpi(mode, 'shepp-logan')
    target = phantom(target_res);
elseif strcmpi(mode, 'flag')
    flag_low = zeros(target_res);
    X = 0:1/(target_res-1):1;
    Y = X;
    rot = 0;
    for i = 1:numel(X)
        flag_low(:,i) = finflag( [ ones(size(X,2),1)*X(i)' Y' ], rot );
    end
    target = flag_low;
elseif strcmpi(mode, 'flag_rot30')
    flag_low = zeros(target_res);
    X = 0:1/(target_res-1):1;
    Y = X;
    rot = 30;
    for i = 1:numel(X)
        flag_low(:,i) = finflag( [ ones(size(X,2),1)*X(i)' Y' ], rot );
    end
    target = flag_low;
end
end