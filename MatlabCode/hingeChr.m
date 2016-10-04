function [ hinges ] = hingeChr( series, span )
%hingeChr identify hinge sites given a vector or a list of vectors
%   series - a matrix containing column vectors
%   span   - smoothing window span applied to `series`

if nargin < 2
    span = 0;
end

[ n, k ] = size(series);

hinges = zeros(k, n);
for i = 1:k
    v = series(:, i)';
    
    if span > 0
        v = smooth(v, span, 'lowess')';
    end
    v = sign(v);

    % obtain hinge sites
    dv = [0 diff(v)];
    dv = abs(sign(dv));
    
    hinges(i, :) = dv;
end

end

