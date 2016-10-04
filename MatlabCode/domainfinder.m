function [ collective_gnm_domains, hinges, gnm_domains_list ] = domainfinder( eigvecs, indices, span )
%DOMAINFINDER identify GNM domains based on given modes
%   eigvecs - eigenvectors used to identify hinge sites
%   indices - indices of nonzero columns/rows
%   span    - smoothing windows span for `eigvecs` to avoid tiny GNM
%             domains caused by noise

factor = 0;

%find domains from each individual mode
[ gnm_domains_list, hinges ]= findGNMDomains(eigvecs, span, 0, factor, indices);
fprintf('finding collective domains...\n')

% find collective domains by taking union of hinge sites
collective_gnm_domains = getCollectiveDomains(gnm_domains_list);

end

