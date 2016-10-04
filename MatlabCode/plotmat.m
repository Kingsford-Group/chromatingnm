function plotmat(mat, varargin)

option = '';
if ~isempty(varargin) 
    option = lower(varargin{1});
end

% if isKirc=1 (input is Kirchhoff matrix), convert all -1s to 0, all 0s to
% 0.5, and all vals greater than 1 to 1 and use imagesc
% if not, just transform all vals to be between 0 and 1 and use imagesc

N = length(mat);
transformedmat = mat;

% Kirchhoff processing
if ~isempty(strfind(option, 'k'))
    transformedmat(mat == 0) = 0.5;
    transformedmat(mat == -1) = 0;
    transformedmat(mat >= 1) = 1;
end

% diagonal mask
di = strfind(option, 'd');
if ~isempty(di)
    di = di(1);
    width = '';
    for i = di + 1: length(option)
        if ~isempty(str2num(option(i)))
            width = strcat(width, option(i));
        end
    end
    if isempty(width); width = '0'; end;
    width = str2num(width);
    I = eye(N);
    for w = 1:width
        diagones = ones(1,N - w);
        I = I + diag(diagones, -w) + diag(diagones, w);
    end
    transformedmat = transformedmat .* (1 - I);
end

% z-score
if ~isempty(strfind(option, 'z'))
    matmean = mean(mat(:));
    matstd  =  std(mat(:));
    transformedmat = (mat - matmean) ./ matstd;
end

% normalization
if ~isempty(strfind(option, 'n'))
    matmax = max(mat(:));
    matmin = min(mat(:));
    transformedmat = (mat - matmin) ./ (matmax - matmin);
end


% masked matrix
if ~isempty(strfind(option, 'm'))
    if length(varargin) > 1 
        maskmat = varargin{2};
        matmax = max(transformedmat(:));
        transformedmat = transformedmat .* (1 - maskmat) + maskmat .* matmax;
    end
end

% plot final matrix
if isempty(strfind(option, '_'))
    figure
end
imagesc(transformedmat)

% truncate caxis
if ~isempty(strfind(option, 't'))
    if length(varargin) >= 2
        factor = varargin{2};
    else
        factor = 3;
    end
    minmat = min(mat(:));
    caxis([minmat abs(minmat) * factor]);
end

% truncate caxis by percentile
if ~isempty(strfind(option, 'p'))
    if length(varargin) >= 2
        factor = varargin{2};
    else
        factor = 1;
    end
    high = prctile(mat(:), 100 - factor);
    low = prctile(mat(:), factor);
    caxis([low high]);
end

cm = 'jet';
if ~isempty(strfind(option, 'c'))
    cm = varargin{2};
end

colormap(cm)
if isempty(strfind(option, '_'))
    colorbar
end