function [ mat, bin, K, fullcov, eigvals0, eigvecs0, cutoff, eigvals, eigvecs, indices ] = readgnm( filename, ZERO )
%READGNM Read GNM results from mat files.
%   mat                - original data
%   bin                - data resolution
%   K                  - Kirchhoff matrix
%   fullcov            - the full covariance matrix
%   eigvals0, eigvecs0 - the raw eigenvalues including zeros
%   eigvals , eigvecs  - the eigenvalues without zeros
%   cutoff             - the used cutoff to construct the Kirchhoff matrix
%   indices            - indices of nonzero rows/columns in raw data

if nargin < 2
    ZERO = 1e-6;
end

% load gnm results
data = load(filename);
mat = data.mat;
bin = data.bin;
K = data.K; 
eigvecs0 = data.eigvecs; 
eigvals0 = data.eigvals; 
fullcov = data.fullcov;
cutoff = data.cutoff;

if isfield(data, 'indices')
    indices = data.indices;
else
    indices = 1:size(eigvecs0, 1);
end

% get meaningful eigenvalues/eigenvectors
i = find(eigvals0 > ZERO, 1);
if ~isempty(i)
    eigvals = eigvals0(i:end);
    eigvecs = eigvecs0(:, i:end);
else
    eigvals = [];
    eigvecs = [];
end
end

