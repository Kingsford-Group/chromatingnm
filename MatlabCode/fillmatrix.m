function [ out ] = fillmatrix( mat, indices )
%FILLMATRIX extend the `mat` with zeros based on indices 
%   This function is used to insert the unmapped regions back to Hi-C map
%   or covariance matrix

matsize = max([max(indices),length(mat)]);
out = zeros(matsize);
out(indices, indices) = mat;
end

