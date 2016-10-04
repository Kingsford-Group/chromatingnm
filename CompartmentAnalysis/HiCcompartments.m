function [compartments, corrmatfull, rawdata] = HiCcompartments(filename)

% input filename should be RaoChr*LAnorm.txt

chrnum = str2double(filename(7));

rawdata = load(filename);

% remove zero rows/columns from rawdata
nonzeros = find(diag(rawdata)~=0);
nonzerodiff = diff(nonzeros);
gapindices = find(nonzerodiff > 1);
gaplist = [nonzeros(gapindices) nonzeros(gapindices+1)];
% zeroindices = find(diag(rawdata)==0);
rawdata_no0 = rawdata(nonzeros,nonzeros);

corrmat = corr(rawdata_no0,rawdata_no0');

%remove diagonal and NaNs
corrmat(1:length(corrmat)+1:length(corrmat(:))) = 0;
corrmat(isnan(corrmat)) = 0;

%use PCA to find boundaries
pcatop5 = HiCpca(corrmat);
if chrnum ~= 4 && chrnum ~=5
    boundarypc = pcatop5(:,1);
else
    boundarypc = pcatop5(:,2);
end

% put back zero rows/columns
firstpcfull = zeros(length(rawdata),1);
firstpcfull(nonzeros) = boundarypc;

corrmatfull = zeros(length(rawdata));
corrmatfull(nonzeros,nonzeros) = corrmat;

%locate compartments by sign changes in first PC
firstpc_zeros = find(firstpcfull(1:end-1).*firstpcfull(2:end) < 0);
compartments = zeros(2*length(firstpc_zeros)+1,2);
compartments(1,:) = [nonzeros(1) firstpc_zeros(1)];
cindex = 2;
for i=2:length(firstpc_zeros)
    gaploc = find(gaplist(:,1) > firstpc_zeros(i-1)+1 & gaplist(:,1) < firstpc_zeros(i));
    if length(gaploc) > 1
        [~,biggestgap] = max(gaplist(gaploc,2)-gaplist(gaploc,1));
        gaploc = gaploc(biggestgap);
    end
    if isempty(gaploc) || gaplist(gaploc,2) - gaplist(gaploc,1) < 10
        compartments(cindex,:) = [firstpc_zeros(i-1)+1 firstpc_zeros(i)];
        cindex = cindex + 1;
    else
        compartments(cindex,:) = [firstpc_zeros(i-1)+1 gaplist(gaploc,1)];
        compartments(cindex+1,:) = [gaplist(gaploc,2)+1 firstpc_zeros(i)];
        cindex = cindex + 2;
    end
end

lastcompartment = find(sum(compartments,2)==0,1);
compartments = compartments(1:lastcompartment,:);
compartments(end,:) = [firstpc_zeros(end) length(corrmatfull)];


% corrmat01 = (corrmat - min(corrmat(:)))/(max(corrmat(:)) - min(corrmat(:)));

% nposcorrs = sum(corrmat(:)>0);
% nnegcorrs = sum(corrmat(:)<0);
% if nposcorrs > nnegcorrs
%     redcols = linspace(.5,1,100);
%     bluecols = linspace(1,.5,floor(nnegcorrs/nposcorrs*100));
% elseif nnegcorrs > nposcorrs
%     bluecols = linspace(1,.5,100);
%     redcols = linspace(.5,1,floor(nposcorrs/nnegcorrs*100));
% else
%     redcols = linspace(.5,1,100);
%     bluecols = linspace(1,.5,100);
% end
% 
% zerovec1 = zeros(length(redcols) + 1,1);
% zerovec2 = zeros(length(bluecols) + 1,1);
% rbcolmap = [[zerovec2; redcols'], zeros(length(zerovec1)+length(zerovec2)-1,1), [bluecols'; zerovec1]];
% 
% redcolmap = ones(101,3);
% redcolmap(2:end,2:3) = 0;
% 
% figure;
% % subplot(121)
% imagesc(log(rawdata))
% colormap(redcolmap)
% title('Normalized Data')
% figure
% % subplot(122)
% imagesc(corrmat)
% colormap(rbcolmap)
% title('Correlation matrix')
