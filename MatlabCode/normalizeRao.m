function [ mat, indices ] = normalizeRao( datafile, normfile, bin, writeTXT )
% handles the normalization of Rao's data

if nargin < 4
    writeTXT = 0;
end

[pathstr,name,ext] = fileparts(datafile);

% load raw data
mat = readRawData(datafile, bin);
% mat = data.mat;
% bin = data.bin;
[n, m] = size(mat);
indices = find(diag(mat) ~= 0);

% read normalization factors
if exist(normfile, 'file')
    [pathstr,name,ext] = fileparts(normfile);
    norm = load(normfile);
    norm = norm(1:n);
    norm_mat = norm * norm';
    mat = mat ./ norm_mat;
    mat(isnan(mat)) = 0;
    mat(isinf(mat)) = 0;
end

% save(fullfile(pathstr, [name, ext, '.mat']), 'mat', 'bin')
if writeTXT
    dlmwrite(fullfile(pathstr, [name, ext, '.txt']), full(mat), '\t')
end
end

