function [ mat ] = readRawData( filename, bin, nrow, ncol)
% read in raw data in sparse format and convert to full matrix

%[pathstr,name,ext] = fileparts(filename);

data = load(filename);
% convert bp to loci
data(:, 1:2) = data(:, 1:2)/bin + 1;

% to guarantee the sparse matrix is a square matrix
if nargin <= 2
    nrow = max(data(:, 1));
    ncol = max(data(:, 2));
    if nrow ~= ncol
        maxidx = max(nrow, ncol);
        data = [data; maxidx maxidx 0];
    end
    % to fill the other half of the matrix
    mat = spconvert(data);
    mat = mat + mat' - diag(diag(mat));
else
    mrow = max(data(:, 1));
    mcol = max(data(:, 2));
    if nrow > mrow || ncol > mcol
        data = [data; nrow ncol 0];
    end
    mat = spconvert(data);
end

% save the matrix to file
% save(fullfile(pathstr, [name, ext, '.mat']), 'mat', 'bin')
end

