function [K, eigvecs0, eigvals0, eigvecs, eigvals, fullcov, indices] = gnmbase(mat, cutoff, option, range)
% calculates Kirchhoff matrix and spectral decomposition, as well as
% nonzero eigenvalues and corresponding eigenvectors

if nargin < 4
    range = 1:size(mat,1);
end

if range == -1
    range = 1:size(mat,1);
end

if nargin < 3
    option = '';
end

mat = mat(range, range);

if strfind(option, 'wozero')
    indices = find(diag(mat) ~= 0);
    mat = mat(indices, indices);
end

% build Kirchhoff matrix with given cutoff
fprintf('building Kirchhoff matrix...\n')
if strfind(option, 'log')
    K = log10(full(mat));
    c = log10(cutoff);
else
    K = full(mat);
    c = cutoff;
end
fprintf('cutoff = %d, option = %s\n', cutoff, option)
K(K <= c) = 0;
K = -K;
for i=1:length(K)
    K(i,i) = 0;
    K(i,i) = -sum(K(i,:));
end


% calculate covariance from first 8 modes
fprintf('decomposing eigenvalues...\n')
[eigvecs0, eigvals0] = eig(K);
eigvals0 = diag(eigvals0);
[eigvals0, I] = sort(eigvals0);
eigvecs0 = eigvecs0(:, I);


if strfind(option, 'lite')
    fullcov = [];
    mat = [];
else
    fprintf('calculating the pseudoinverse...\n')
    fullcov = pinv(K);
end


ZERO = 1e-6;
% get meaningful eigenvalues/eigenvectors
i = find(eigvals0 > ZERO, 1);
if ~isempty(i)
    eigvals = eigvals0(i:end);
    eigvecs = eigvecs0(:, i:end);
else
    eigvals = [];
    eigvecs = [];
end
fprintf('Done.\n')
end