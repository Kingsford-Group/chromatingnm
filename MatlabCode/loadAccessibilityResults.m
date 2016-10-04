function [cum, mean, count, startIDX, endIDX, data ] = loadAccessibilityResults( filename, bin )
%loadAccessibilityResults loads the ATAC or DNase data file in mat format

data = load(filename);
starts = data.starts;
ends = data.ends;
signal = data.signal;

if nargin < 2
    bin = 40000;  % 40kb
end

startIDX = ceil((starts + 1)/bin);
endIDX = ceil(ends/bin);

series = zeros(1, max(endIDX));
count = zeros(1, max(endIDX));

for i = 1:length(startIDX)
    series(startIDX(i)) = series(startIDX(i)) + signal(i);
    count(startIDX(i)) = count(startIDX(i)) + 1;
end
cum = series;
mean = series ./ count;
mean(isnan(mean)) = 0;
end

