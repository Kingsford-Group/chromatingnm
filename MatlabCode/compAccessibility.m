function [correlations, msfvec, dnasevec, atacvec] = compAccessibility( GNM_file, ATAC_file, DNase_file, n_modes, span)
%compATAC compares GNM to open chromatin data.
%   filename: the filename of GNM results 
%   ATAC_file, DNase_file: open chromatin data file
%   n_modes: the number of modes used to calculate the square fluctuations


% load gnm results
vars = whos('-file', GNM_file);
has_sqf = ismember('sqf', {vars.name});
if has_sqf
    [ sqf, modes, bin, cutoff ] = readsqf( GNM_file );
    maxmodes = length(modes);
    indices = 1:length(sqf);
else
    [ ~, bin, ~, fullcov, ~, ~, cutoff, eigvals, eigvecs, indices ] = readgnm( GNM_file );
    maxmodes = length(eigvals);
end
%[eigvecs, eigvals] = modefilter(eigvecs, eigvals);
xlbl = sprintf('loci (%dkb)', bin/1e3);

% load ATAC results
if nargin > 1
    [ATACs, ~, ~, ATACStarts, ATACEnds, data] = loadAccessibilityResults(ATAC_file, bin);
else
    ATACStarts = [];
    ATACEnds = [];
end

% load DNaseI results
if nargin > 2
    [DNase, ~, ~, DNaseStarts, DNaseEnds, data] = loadAccessibilityResults(DNase_file, bin);
else
    DNaseStarts = [];
    DNaseEnds = [];
end

if nargin < 4
    n_modes = 10;
end

if nargin < 5
    span = 200e3/bin + 1;
end

if n_modes < 0 || n_modes > maxmodes
    n_modes = maxmodes;
end

% calculate covariance from first n modes without screening
if ~has_sqf
    [sqf, ~] = calcSqfluct(eigvecs, eigvals, 1:n_modes);
    sqf = prepseries(sqf, 'fill', indices);
%     fprintf('MSFs from first %d modes\n', n_modes)
end

% preprocessing the series
msfvec = prepseries(sqf, 'norm');
dnasevec = prepseries(DNase', 'norm', 'smooth', span);
atacvec = prepseries(ATACs', 'norm', 'smooth', span);
exps = {'ATAC', 'DNase'};


% trim the series into the same length
mutlen = min([max(indices), length(dnasevec), length(atacvec)]);
maxlen = max([length(msfvec), length(dnasevec), length(atacvec)]);

corr1 = corr(msfvec(1:mutlen), dnasevec(1:mutlen), 'type', 'Spearman');
corr2 = corr(msfvec(1:mutlen), atacvec(1:mutlen), 'type', 'Spearman');
corr3 = corr(dnasevec(1:mutlen), atacvec(1:mutlen), 'type', 'Spearman');
correlations = [corr1 corr2 corr3];
% fprintf('spearman correlation of sqfluct and %s: %6.4f\n', exps{1}, corr1)
% fprintf('spearman correlation of sqfluct and %s: %6.4f\n', exps{2}, corr2)
% fprintf('spearman correlation of ATAC and DNase: %6.4f\n', corr3)

