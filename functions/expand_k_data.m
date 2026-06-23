function k_data_ext = expand_k_data(k_data_primitive, mode, param)
% EXPAND_K_DATA Generate extended frequency set from primitive vectors.
%
% Each primitive vector k = (k1,k2) is expanded by integer multiples m*k.
% Only multiples that remain inside the specified domain are kept:
%   - 'box'  : |k1|, |k2| <= N2
%   - 'ball' : ||k|| <= R
% The associated data Dv is copied for all valid multiples.
%%
k_data_ext = [];

for i = 1:size(k_data_primitive,1)
    
    v = k_data_primitive(i,1:2);
    
    % Skip (0,0)
    if v(1) == 0 && v(2) == 0
        continue
    end
    
    Dv = k_data_primitive(i,3:end);
    
    % How many multiples maximum inside the bigger ball or rectangle
    if strcmp(mode,'box')
        N2 = param;
        m_max = floor(N2 / max(abs(v)));
    elseif strcmp(mode,'ball')
        R = param;
        m_max = floor(R / norm(v));
    end
    
    % Loop through multiples
    for m = 1:m_max
        
        k = m * v;
        
        % Restriction
        if strcmp(mode,'box')
            if abs(k(1)) > N2 || abs(k(2)) > N2
                continue
            end
        elseif strcmp(mode,'ball')
            if norm(k) > R
                continue
            end
        end
        
        k_data_ext = [k_data_ext; k Dv];
    end
end

end
