function [ccdds, medianGNMweights] = findCCDDs(fullcov,outputname,nonzeroindices)

if nargin < 2 || isempty(outputname)
    nofileoutput = 1;
else
    nofileoutput = 0;
end

if nargin < 3 || isempty(nonzeroindices) || length(nonzeroindices) < length(fullcov)
    nonzeroindices = 1:length(fullcov);
end

res = 5000;

% binarize matrix with threshold of |min(covmatrix)|
theta = abs(min(fullcov(:)));
fullcovbw = fullcov;
fullcovbw(fullcov > theta) = 1;
fullcovbw(fullcov < theta) = 0;
fullcovbw = tril(fullcovbw);

connstruct = bwconncomp(fullcovbw);
numPixels = cellfun(@numel,connstruct.PixelIdxList);

mediancovval = @(pixlist) median(fullcov(pixlist));
medianofRegion = cellfun(mediancovval, connstruct.PixelIdxList);

[numPixels_sorted, sortidxnum] = sort(numPixels,'descend');

[diag_ind_x, diag_ind_y] = ind2sub(size(fullcov),connstruct.PixelIdxList{sortidxnum(1)});

% find widest point along diagonal for distance threshold
distthresh = max(abs(diag_ind_x - diag_ind_y));

stopidx = find(numPixels_sorted < 2,1)-1; % ignore singletons

% fprintf(strcat('Total number of regions: ',num2str(stopidx-1),'\n'))

fullcovbw = fullcovbw + fullcovbw';
ccdds = zeros(stopidx,4);
medianGNMweights = zeros(stopidx,1);
ndoms = 0;
for i=2:stopidx
    [x, y] = ind2sub(size(fullcov),connstruct.PixelIdxList{sortidxnum(i)});
    if abs(y-x) > distthresh %check that CCDD is outside main diagonal
        ndoms = ndoms + 1;
%         fprintf(strcat('Finding boundaries of region ',num2str(ndoms),'\n'))

        maxrectarea = 0;
        % find maximum area rectangle contained within region
        for j=1:length(x)
            if length(x)-j+1 > maxrectarea
                brcorner = find(fullcovbw(y(j):end,x(j)) == 0,1);
                if isempty(brcorner)
                    brcorner = length(fullcovbw) - y(j) + 1;
                end
                tlcorner = find(fullcovbw(y(j),x(j):end) == 0,1);
                if isempty(tlcorner)
                    tlcorner = length(fullcovbw) - x(j) + 1;
                end
                rectarea = (brcorner-1)*(tlcorner-1);
                rect = fullcovbw(y(j):y(j)+brcorner-2,x(j):x(j)+tlcorner-2);
                if ~isempty(find(rect == 0,1));
                    continue
                elseif rectarea > maxrectarea
                    maxrectarea = rectarea;
                    biggestrect = [y(j), x(j); y(j)+brcorner-1, x(j)+tlcorner-1; rectarea, 0];
                end
            end
        end
        
        ymin = nonzeroindices(biggestrect(1,1));
        xmin = nonzeroindices(biggestrect(1,2));
        ymax = nonzeroindices(biggestrect(2,1));
        xmax = nonzeroindices(biggestrect(2,2));
        
        ccdds(ndoms,:) = [(xmin-1)*res, xmax*res-1, (ymin-1)*res, ymax*res-1];
        medianGNMweights(ndoms) = medianofRegion(sortidxnum(i));
    end
end

ccdds(ndoms+1:end,:) = [];
medianGNMweights(ndoms+1:end) = [];

if ~nofileoutput
    filename = outputname;
    dlmwrite(filename, [ccdds, medianGNMweights],'delimiter','\t','precision',10);
end

