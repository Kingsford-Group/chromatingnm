function [ gnm_domains_list, hinges ] = findGNMDomains( series, smoothingspan, min_domain_length, min_domain_mean_factor, indices)
%findGNMDomains
%   series can be a matrix where each column is a eigenvector.
%   E.g. [ gnm_domains_list, hinges ]= findGNMDomains(eigvecs(:,1:n_modes), span, 0, factor, indices);

if nargin < 3
    min_domain_length = 2;
end

if nargin < 4
    min_domain_mean_factor = 0.05;
end

if nargin < 5
    indices = 1:size(series, 1);
end

gnm_domains_list = {};

% identify breaks
hinges = hingeChr(series, smoothingspan);
hinges = prepseries(hinges', 'fill', indices)';
series = prepseries(series, 'fill', indices);
hinges2 = hinges;
hinges2(:, 1) = 1;
hinges2(:, end) = 1;
for i = 1:size(hinges2, 1)
    sites = find(hinges2(i, :) > 0);
    gnm_domains_raw = zeros(length(sites)-1, 2);
    % construct consecutive domains from hinge sites 
    for j = 1:length(sites)-1
        gnm_domains_raw(j, :) = [sites(j) sites(j+1)];
    end
    
    msfs = [];
    for j = 1:size(gnm_domains_raw, 1)
        domain = gnm_domains_raw(j, :);
        sf = mean(abs(series(domain(1):domain(2))));
        msfs = [msfs sf]; %#ok<*AGROW>
    end

    gnm_domains = [];
    for j = 1:size(gnm_domains_raw, 1)
        domain = gnm_domains_raw(j, :);
        length_domain = domain(2) - domain(1) + 1;
        sf = msfs(j);
        if length_domain >= min_domain_length && sf >= std(series(:, i)) * min_domain_mean_factor
            gnm_domains =[gnm_domains; domain];
        end
    end
    
    gnm_domains_list{i} = gnm_domains;
end

end

