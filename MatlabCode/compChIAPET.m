function [ ChIA_cov, non_ChIA_cov ] = compChIAPET( GNM_file, ChIA_file, chr, min_dist, sup_plot)


%close all
if nargin < 4
    min_dist = 0;
end

if nargin < 5
    sup_plot = 0;
end
[ mat, bin, ~, fullcov, ~, ~, cutoff, eigvals, eigvecs, indices ] = readgnm( GNM_file );
[pos1, pos2, chrom1, chrom2, pair_name, Z, ~, ~, ~, ~, ~] = readChIAPET(ChIA_file, bin, chr);

%% plot covariance matrix with locations of ChIA-pet and background pairs
cov = fillmatrix(fullcov, indices);
mincov = min(cov(:));

if ~sup_plot
    covwodiag = cov - diag(diag(cov));
    figure()
    set(gcf,'Renderer','OpenGL');
    plotmat(covwodiag, 'c_', 'jet');
    caxis([mincov abs(mincov) * 5])
    % caxis([mincov maxcov])
    colorbar()
end

%% obtain ChIA covariance
ChIA_cov = [];
ChIA_pos = [];
ChIA_mat = [];

for i = 1:length(pos2)
    start1 = pos1(i,1);
    end1 = pos1(i,2);
    start2 = pos2(i,1);
    end2 = pos2(i,2);
    for p1 = start1:end1
        for p2 = start2:end2
            if abs(p1 - p2) >= min_dist
                ChIA_cov = [ChIA_cov cov(p1, p2)]; %#ok<*AGROW>
                ChIA_pos = [ChIA_pos; p1, p2]; 
                ChIA_mat = [ChIA_mat mat(p1, p2)]; 
            end
        end
    end
end
% make sure that pos1 is greater than pos2 for convenience
ChIA_pos = sort(ChIA_pos, 2, 'descend');

if ~sup_plot
    hold on
    scatter(ChIA_pos(:,1), ChIA_pos(:,2), 'w.')
    hold off
end

%% sample background covariance
non_ChIA_pos = [];
non_ChIA_cov = [];
non_ChIA_mat = [];


for i = 1:size(ChIA_pos, 1)
    p1 = ChIA_pos(i, 1); p2 = ChIA_pos(i, 2);
    % obtain the distance of the ChIAs
    d = p1 - p2;
    
    % obtain the position with the same distance but on the other side
    % p4-----p2-----p1-----p3
    p3 = p1 + d;
    lia = ismember([p1 p3; p3 p1], ChIA_pos, 'rows');
    if ~any(lia) && p3 <= length(cov)
        non_ChIA_pos = [non_ChIA_pos; p1 p3];
        non_ChIA_cov = [non_ChIA_cov cov(p1, p3)];
        non_ChIA_mat = [non_ChIA_mat mat(p1, p3)];
    end
    
    % obtain the position with the same distance but on the other side
    p4 = p2 - d;
    lia = ismember([p2 p4; p4 p2], ChIA_pos, 'rows');
    if ~any(lia) && p4 > 0
        non_ChIA_pos = [non_ChIA_pos; p4 p2];
        non_ChIA_cov = [non_ChIA_cov cov(p2, p4)];
        non_ChIA_mat = [non_ChIA_mat mat(p2, p4)];
    end
end

if ~sup_plot
    hold on
    scatter(non_ChIA_pos(:,1), non_ChIA_pos(:,2), 'r.')
    hold off
end

%% plot the histogram comparing ChIA covariance to background covariance

mincov = min([non_ChIA_cov ChIA_cov]);
maxcov = max([non_ChIA_cov ChIA_cov]);
bins = mincov:3e-7:maxcov;

if ~sup_plot
figure
hold on
histogram(non_ChIA_cov, bins, 'normalization', 'probability');
histogram(ChIA_cov, bins, 'normalization', 'probability');
legend('background', 'ChIA')
xlabel('Covariance')
ylabel('Probability')
mean_ChIA = mean(ChIA_cov);
mean_non_ChIA = mean(non_ChIA_cov);
plot([mean_non_ChIA mean_non_ChIA], ylim(), 'b--');
plot([mean_ChIA mean_ChIA], ylim(), 'r--');

[h, p] = ttest2(non_ChIA_cov, ChIA_cov, 'Vartype','unequal');
fprintf('interaction frequency h = %d, p = %d\n', h, p)
hold off
end

end

