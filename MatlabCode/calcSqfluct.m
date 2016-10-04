function [ sqf, outmodes ] = calcSqfluct( eigvecs, eigvals, modes )
%calcSqfluct calculates the MSFs from eigenvalues and vectors based on a subset of modes.
% eigvecs & eigvals - full set of eigenpairs obtained by GNM
% modes             - a vector of mode indices

[n, m] = size(eigvecs);
if nargin < 3
    modes = 1:m;
end


% calculate covariance from first n modes
sqf = zeros(n, 1);
outmodes = 0;

j = 0;
for i = modes
    outmodes = outmodes + 1;
    sqfi = (eigvecs(:,i).*eigvecs(:,i))/eigvals(i);
    sqf = sqf + sqfi;
    
    j = j + 1;
end

end

