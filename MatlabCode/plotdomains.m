function [ p ] = plotdomains( domains_list, upper, varargin )
%PLOTDOMAINS 

if nargin < 2
    upper = 1;
end

if nargin < 3
    varargin = {};
end
hold on
for item = domains_list
    domains = item{1};
    for i = 1:size(domains, 1)
        domain = domains(i, :);
        if upper
            p = plot([domain(1) domain(2) domain(2)], [domain(1) domain(1) domain(2)], varargin{:});
        else
            p = plot([domain(1) domain(1) domain(2)], [domain(1) domain(2) domain(2)], varargin{:});
        end
    end
end
hold off
end

