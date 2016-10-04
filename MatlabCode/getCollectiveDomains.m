function [ collective_domains ] = getCollectiveDomains( domain_cells, min_domain_length )
%getCollectiveDomains 
%   domain_cells contains a list of domains identified by different
%   methods/resolutions.

if nargin < 2
    min_domain_length = 1;
end

% Preprocessing. Restore the domain starts and ends into seperate lists.
% Find the max ends of domains.
n_domain_cells = length(domain_cells);
domain_starts = [];
domain_ends = [];
series_lengths = [];
for i = 1:n_domain_cells
    domains = domain_cells{i};
    series_lengths = [series_lengths; domains(end)]; %#ok<*AGROW>
    domain_starts = [domain_starts; domains(:, 1)];
    domain_ends = [domain_ends; domains(:, 2)];
end

series_length = max(series_lengths);

collective_domains = [];
curr_start = nan;
for i = 1:series_length
    is_start = ~isempty(find(domain_starts == i, 1));
    is_end = ~isempty(find(domain_ends == i, 1));
    
    if is_start % if current position is labelled as a START in any domains
        if isnan(curr_start) % if last domain is finished, then create a new domain
            
        else  % if last domain is not finished, then finish it and create a new domain
            collective_domains = [collective_domains; curr_start, i-1];
        end
        curr_start = i;
    elseif is_end % if current position is labelled as an END in any domains
        if isnan(curr_start) % if last domain is finished
            collective_domains = [collective_domains; i-1, i];
        else  % if last domain is not finished, then finish it
            collective_domains = [collective_domains; curr_start, i];
        end
        
        % check if current location is within a domain
        in_domain = 0;
        for j = 1:length(domain_starts)
            startj = domain_starts(j);
            endj   = domain_ends(j);
            if i > startj && i < endj 
                in_domain = 1;  % the current location is found to be in a domain
                break           % then escape the checking procedure
            end
        end
        
        if in_domain          % if in the domain, 
            curr_start = i;   % then create a new domain after terminating from current domain
        else                  % if not, then simply terminate current domain
            curr_start = nan;
        end
    end
end

temp_collective_domains = collective_domains;
collective_domains = [];
for i = 1:size(temp_collective_domains, 1)
    domain = temp_collective_domains(i, :);
    length_domain = domain(2) - domain(1) + 1;
    if length_domain >= min_domain_length 
        collective_domains =[collective_domains; domain];
    end
end
end

