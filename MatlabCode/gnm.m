function [eigvecs, eigvals, fullcov, indices, gnmfile] = gnm(datafile, normfile, bin, cutoff, option, outputpath)
% GNM performs gnm analysis given 2d contact map
% datafile (string): the raw contact map file.
% normfile (string): the normalization file.
% bin (int): bin size, e.g. 5e3 for 5kb resolution.
% cutoff (int): the cutoff used in constructing Kirchhoff matrix.
% option (string): "wozero" removes the unmapped regions;
%
%                  "log" for log transformation of the data;
%
%                  "lite" for smaller output file. With this option on, the
%                  function will not write the original matrix into the
%                  output file;
%
%                  "writetxt" is to also output the matrix in txt format.
%                  Note that it will have the unmapped regions even if you
%                  have the "wozero" option;
%
%                  The options can be used sequentially, e.g.
%                  "wozero.writetxt" will output the txt file and remove
%                  the unmapped region.

close all

if nargin < 3
    option = '';
end

if ~isempty(option) && ~strcmp(option(1), '.') 
    option = ['.' option];
end

if strfind(option, 'writetxt')
    writeTXT = 1;
    option = strrep(option, '.writetxt', '');
else
    writeTXT = 0;
end

fprintf('preparing the raw data...\n')
[mat, ~] = normalizeRao(datafile, normfile, bin, writeTXT);
[K, eigvecs, eigvals, ~, ~, fullcov, indices] = gnmbase(mat, cutoff, option);


[pathstr,name,ext] = fileparts(datafile);
if exist(normfile, 'file')
    [pathstr,name,ext] = fileparts(normfile);
end

if nargin == 6
    pathstr = outputpath;
end

output = sprintf('%s.c=%d%s', [name, ext], cutoff, option);
output = strcat(output, '.gnm', '.mat');
fprintf('saving results...\n')
gnmfile = fullfile(pathstr, output);
save(gnmfile, 'mat', 'bin', 'K', 'eigvecs', 'eigvals', 'fullcov', 'cutoff', 'indices', '-v7.3')
fprintf('Done.\n')
end