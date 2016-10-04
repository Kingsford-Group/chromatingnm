function [ cov ] = covn( eigvecs, eigvals, n_modes )
%calculate covariance matrix for given number of modes or subset of eigenpairs
% eigvecs & eigvals - full set of eigenpairs obtained by GNM
% n_modes           - a vector of mode indices or a scalar indicating number of modes
%                     that should be used

[n, m] = size(eigvecs);
if nargin < 3
    n_modes = m;
elseif n_modes < 0
    n_modes = m;
end

% calculate covariance from first n modes
cov = zeros(n, n);

if length(n_modes) > 1
    mode_range = n_modes;
else
    mode_range = 1:n_modes;
end

for i = mode_range
    modecontribution = (eigvecs(:,i)*eigvecs(:,i)')/eigvals(i);
    cov = cov + modecontribution;
end

end

