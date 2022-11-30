%load('26jul22workspace.mat')
% cluster on whole group 
%it's got everybody, clinical groups and hetero group
% do PCA
% k: 2-10
%then calculate ARIs between these groups and the clinically homogeneous
%groups (and then again with the clusters)
        [~, score, ~, ~, explained]= pca(Zs);
        esum = zeros(88,1);
        esum(1) = explained(1);
        for b = 2:88
            esum(b) = esum(b-1) + explained(b);
        end
        f = find(esum > 75);
        fs = f(1,1);
        [~, ka, kms] = clusk(score(:,1:f(1,1)), 10); 
        
        % maybe I want to save score, explained, and kms? and ka?
%% ARI clinical groups
for k = 1:10
    ri_clingroups{k} = rand_index(kms(:,k), G, 'adjusted');
end

%% clusters
itotal = {ia, im, is, ic, il, ise};
for i = 1:6
    for k = 1:10
        ri_clus{i,k} = rand_index(kms(itotal{i},k), kms_final{i}, 'adjusted');
    end
end
%% next
figure
bar(explained(1:32,:))
xlabel('PCA Component')
ylabel('% of Total Variance Explained')
legend('Anhedonia', 'Depressed Mood', 'Somatic Disturbance', 'Chronicity', 'Late Onset', 'Acute Impairment', 'Heterogeneous')
set(gca, 'FontSize', 14)





%% Functions
function [idx_kms, a, kms] = clusk(Z, varargin)
for n = 1:varargin{1}
    kms(:,n) = kmeans(Z, n);
    kms_sil(:,n) = mean(silhouette(Z, kms(:,n)));
end
[~, a] = max(kms_sil);
idx_kms = kms(:,a);
end