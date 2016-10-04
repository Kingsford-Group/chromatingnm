%% Script to run full GNM analysis
% This is a script that reproduces work presented in ...

%% Prepare files and calculate GNM components
chrstr = '22';

rawhicfile = strcat('InputData/chr',chrstr,'_5kb.RAWobserved');
normfile = strsplit(rawhicfile,'.');
normfile = strcat(normfile{1},'.VCnorm');
res = 5e3;

% first calculate results of spectral decomposition
[eigvecs, eigvals, fullcov_nozeros, nonzeroindices, gnmfile] = gnm(rawhicfile, normfile, res, 0, 'wozero', 'OutputData/');

% for a covariance matrix from a subset of modes
n_modes = 500;
cov_nmodes = covn(eigvecs,eigvals,n_modes);

%% Comparison with DNase/ATAC-seq
dnasefile_raw = 'InputData/GM12878DukeDNaseSeq.pk';
atacfile_raw = 'InputData/GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed';
% read .pk and .bed files and write to .mat files
readNarrowPk(dnasefile_raw, 1, 'OutputData/');
readBED(atacfile_raw, 1, 'OutputData/');

[~,fileprefix,~] = fileparts(dnasefile_raw);
dnasefile_mat = strcat('OutputData/',fileprefix,'.chr',chrstr,'.mat');
[~,fileprefix,~] = fileparts(atacfile_raw);
atacfile_mat = strcat('OutputData/',fileprefix,'.chr',chrstr,'.mat');

span = 2e5/res;
[correlations, msfvec, dnasevec, atacvec] = compAccessibility(gnmfile, atacfile_mat, dnasefile_mat, n_modes, span);
%correlations vector: Spearman correlation of [ MSFwithDNase, MSFwithATAC, DNasewithATAC]

%% Identify GNM domains
eigvec_span = 15;
n_modes_for_domains = 10;
[ collective_gnm_domains, hinges, gnm_domains_list ] = domainfinder(eigvecs(:,1:n_modes_for_domains), nonzeroindices, eigvec_span);

% plot GNM domains on top of covariance matrix
fullcov = fillmatrix(fullcov_nozeros, nonzeroindices);
plotmat(fullcov, 'p') % plot the full covariance matrix
hold on
plotdomains({collective_gnm_domains}, 0, 'y:', 'LineWidth', 2); % the upper triangle of the domain indicators
plotdomains({collective_gnm_domains}, 1, 'y:', 'LineWidth', 2); % the lower triangle
hold off
title('GNM domains')

%% Compare GNM domains with TADs and compartments

% load in compartment boundaries
compartments = load('InputData/RaoChr22_5kb_compartments.txt');
% details on calculating compartments can be found in Lieberman-Aiden et
% al., 2009

% load in TAD boundaries, calculated by Armatus
tadfileid = fopen('InputData/RaoChr22KR_5kb.consensus.txt');
taddata = textscan(tadfileid, '%s %d %d');
taddata = taddata(2:3);
tads = [taddata{1}, taddata{2}]./res;
fclose(tadfileid);

% calculate VI between compartments, TADs, and GNM domains
comp_gnmdom_vi = vi(compartments, collective_gnm_domains, length(fullcov));
tad_gnmdom_vi = vi(tads, collective_gnm_domains, length(fullcov));

% visualization of domain comparisons
plotmat(fullcov, 'p') % plot the full covariance matrix
hold on
plotdomains({collective_gnm_domains}, 0, 'y:', 'LineWidth', 2); % the upper triangle of the domain indicators
plotdomains({compartments}, 1, 'r:', 'LineWidth', 2); % the lower triangle
legend('GNM domains','Compartments')
hold off

plotmat(fullcov, 'p') % plot the full covariance matrix
hold on
plotdomains({collective_gnm_domains}, 0, 'y:', 'LineWidth', 2); % the upper triangle of the domain indicators
plotdomains({tads}, 1, 'r:', 'LineWidth', 2); % the lower triangle
legend('GNM domains','TADs') %fix legends
hold off

%% Validation through ChIA-PET pairs

[ChIA_pet_covvalues, background_covvalues] = compChIAPET(gnmfile, 'InputData/ENCFF002EMO.tsv',22,0,1);

%% Cross-correlated distal domains (CCDDs)

outputfilename = strcat('OutputData/ccddLocs_chr',chrstr,'.txt');
[ccdds, medianGNMweights] = findCCDDs(fullcov_nozeros, outputfilename,nonzeroindices);

% co-expression analysis done outside of Matlab, see README for details.