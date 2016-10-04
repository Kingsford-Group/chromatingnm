function [ out ] = prepseries( in, varargin )
%PREPSERIES pre-process input data series
%   in - a matrix of data series. Each column is one series.

function [ out ] = prep1series( in, varargin )
out = in;
skip = 0;
if ~isempty(varargin) 
    for i =1:length(varargin)
        if ~skip
            arg = varargin{i};
            if strfind(arg, 'norm')
                out = out/norm(out);
            elseif strfind(arg, 'smooth')
                out = smooth(out, varargin{i+1});
                skip = 1;
            elseif strfind(arg, 'fill')
                indices = varargin{i+1};
                temp = zeros(max(indices), 1);
                temp(indices) = out;
                out = temp;
                skip = 1;
            end
        else
            skip = 0;
        end
    end
end
end

n_cols = size(in, 2);
out1 = prep1series(in(:, 1), varargin{:});
out  = zeros(length(out1), n_cols);
out(:, 1) = out1;
for c = 2:n_cols
    out(:, c) = prep1series(in(:, c), varargin{:});
end
end

