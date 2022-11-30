%%PCA of Zs

% add paths
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data/allage')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Code')
%load Z_TNM_case.csv
Zs = table2array(readtable('Z_TNM_case.csv'));
%load Z_TNM_case_eid.csv
Zeids = table2array(readtable('Z_TNM_case_eid.csv'));

%%load clinical group eids
anhedonia_eid = table2array(readtable('HomogeneousGroups_anhedonia.csv'));
mood_eid = table2array(readtable('HomogeneousGroups_mood.csv'));
somatic_eid = table2array(readtable('HomogeneousGroups_somatic.csv'));
chronic_eid = table2array(readtable('HomogeneousGroups_chronic.csv'));
lateonset_eid = table2array(readtable('HomogeneousGroups_lateonset.csv'));
severe_eid = table2array(readtable('HomogeneousGroups_severe.csv'));
hetero_eid = table2array(readtable('HomogeneousGroups_HeterogeneousComparison.csv'));
% anhdep_eid = table2array(readtable('HomogeneousGroups_anhdep.csv'));

%%select Zs for groups
[~, ia] = intersect(Zeids(:,1), anhedonia_eid(:,1)); Z{1} = Zs(ia, :);
[~, im] = intersect(Zeids(:,1), mood_eid(:,1)); Z{2} = Zs(im, :);
[~, is] = intersect(Zeids(:,1), somatic_eid(:,1)); Z{3} = Zs(is, :);
[~, ic] = intersect(Zeids(:,1), chronic_eid(:,1)); Z{4} = Zs(ic, :);
[~, il] = intersect(Zeids(:,1), lateonset_eid(:,1)); Z{5} = Zs(il, :);
[~, ise] = intersect(Zeids(:,1), severe_eid(:,1)); Z{6} = Zs(ise, :);
[~, ih] = intersect(Zeids(:,1), hetero_eid(:,1)); Z{7} = Zs(ih, :);
% [~, iad] = intersect(Zeids(:,1), anhdep_eid(:,1)); Zc{1} = Zs(iad, :);

%% rerun ANOVA
%%%create G 
% [~, ih] = intersect(Zeids(:,1), hetero_eid(:,1));
G(ia, 1) = 1;
G(im, 1) = 2;
G(is, 1) = 3;
G(ic, 1) = 4;
G(il, 1) = 5;
G(ise, 1) = 6;
G(ih, 1) = 7;
for n = 1:size(Zs,2)
    [anovaPs(n), tbl{n}, stats{n}] = anovan(Zs(:,n),G,'display','off');
%     c{n} = multcompare(stats);
end
[pthr,pcor,padj] = fdr(anovaPs);
sigidx = find(pcor < 0.0500);
for i = 1:length(sigidx)
c{i} = multcompare(stats{sigidx(i)});
h{i}= find(c{i}(:,6) < 0.0500);
end

%clear ia im is ic il ise anhedonia_eid mood_eid somatic_eid chronic_eid lateonset_eid severe_eid Zs Zeids

%% bootstrap ANOVA
%randomly select with replacement
% run ANOVA
%save what's significant (not post hoc though for now)
for i = 1:100
    ix{i} = randi(length(Zs), [length(Zs),1]);
    data = Zs(ix{i},:);
    group = G(ix{i},:);
for n = 1:size(Zs,2)
    anovaPs(n) = anovan(data(:,n),group,'display','off');
end
[pthr,pcor,padj] = fdr(anovaPs);
sigidx100{i} = find(pcor < 0.0500);

end

%% visualize
% color = repmat([0.3010 0.7450 0.9330], [88,1]);
% for i = 1:16
% color(sigidx(i),:) = [0.6350 0.0780 0.1840];
% end
% histogram(cell2mat(sigidx100), 88, 'FaceColor', color)
% ylabel('Number of Bootstraps Imaging Feature Was Significant')
% xlabel('Imaging Feature (1-88)')

figure;
h = histogram(cell2mat(sigidx100), 88);
b = bar(1:88,h.Values);
b.FaceColor = 'flat';
for i = 1:16
b.CData(sigidx(i),:) = [.5 0 .5];
end
ylabel('Number of Bootstraps Imaging Feature Was Significant')
set(gca,'FontSize',12);
% xlabel('Imaging Feature (1-88)  **true significant imaging features in purple')
xticks(1:88)
xticklabels(IDPsselect89order1)
ax = gca;
ax.XAxis.FontSize = 9;
set(gca,'view',[90 -90])
%% Intraclass correlation
for i = 1:7
[r{i}, LB{i}, UB{i}, F{i}, df1{i}, df2{i}, p{i}] = ICC(Z{i}, 'C-k'); %https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc
end

%% kms bootstrap on PCA
for i = 1:1000
    for n = 1:6
        ix = randi(length(Z{n}), [length(Z{n}),1]);
        data = Z{n}(ix,:);
        [~, score, ~, ~, explained]= pca(Z{n});
        esum = zeros(88,1);
        esum(1) = explained(1);
        for b = 2:88
            esum(b) = esum(b-1) + explained(b);
        end
        f = find(esum > 75);
        fs(i) = f(1,1);
        [kms{i,n}, ka(i,n)] = clusk(score(:,1:f(1,1)), 10); 
        i
    end
end
%%find optimal cluster number
% for i = 1:6
%     cluster(i) = mode(ka(:,i));
% end
%%
%%ARIscore_final
folds = 10;
for x = 1:100
    for i = 1:folds
        for n = 1:6
            c = cvpartition(size(Z{n},1),'KFold',folds);
            idxTrain{x,i,n} = training(c, i);
            Z_Train{x,i, n} = Z{n}(idxTrain{x,i,n},:);
            idxNew{x,i,n} = test(c, i);
            Z_New{x,i, n} = Z{n}(idxNew{x,i,n},:);
        end
    end

    for i = 1:folds
        for n = 1:6
            [coeff{x,i,n}, score_temp{x,i,n}, ~,~,explained] = pca(Z_Train{x, i, n}); %%I think I need to save out the weights of this PCA to use 2 steps down
            esum = zeros(88,1);
            esum(1) = explained(1);
            for b = 2:88
                esum(b) = esum(b-1) + explained(b);
            end
            f = find(esum > 75);
            %%%%I don't understand why centroid is 2x27 when score_temp is
            %%%%282x32
            [kms_temp{x,i,n}, centroid{x,i,n}] = kmeans(score_temp{x,i,n}(:,1:f(1,1)), cluster(n)); 
            score = Z_New{x,i, n}*coeff{x,i,n}; 
            score_new{x,i,n} = score(:,1:f(1,1));
        end
    end
    
    for i = 1:folds
        for n=1:6
            for j = 1:size(centroid{x,i,n}, 1)
                for k = 1:size(score_new{x,i,n},1)
                    A{x,i,n}{j,k} = sum((centroid{x,i,n}(j,:)-score_new{x,i,n}(k,:)).^2,2); %j = number of clusters, k = number of subjects
                    [~,cl_idx{x,i,n}(k)] = min(cell2mat(A{x,i,n}(:,k)));
                end
            end
            total{x,i, n} = zeros(1,length(idxTrain{x,i,n}));
            total{x,i, n}(:,idxTrain{x,i,n}) = kms_temp{x,i,n};
            total{x,i, n}(:,idxNew{x,i,n}) = cl_idx{x,i,n};     
        end
    end
end
%%
     %%%ARI
 for n = 1:6
    for x = 1:100
        for i = 1:folds
            for j = i+1:folds
                ARI_val(i,j) = rand_index(total{x,i,n},total{x,j,n}, 'adjusted'); %https://www.mathworks.com/matlabcentral/fileexchange/49908-adjusted-rand-index
                ARI_all{n,x}(i,j) = ARI_val(i,j);

            end  
        end
        ARI_val100 = mean(ARI_val(find(~tril(ones(size(ARI_val))))));
    end
    ARI_dist{n} = mean(ARI_val100);
    n
 end

%% save out important variables
for n = 1:6
     [~, score_final{n}]= pca(Z{n});
     kms_final{n} = kmeans(score_final{n}, cluster(n)); 
%    kms_final{n} = kmeans(score_final{n}(:, 1:choice), cluster(n)); 
end
save('PCA_ARI_allage_75jul8.mat', 'ARI_dist', 'total', 'kms', 'ka', 'cluster')
save('finalPCAscores_kms_allage_75jul8.mat', 'kms_final', 'score_final')


%% validations

%% cognition
%load cognition again. select the 2, multiply the 2 by the weights
D = readtable('Cognition.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
Reasoning = [D.x20016_0_0 D.x20016_2_0];
ReactionTime = [D.x20023_0_0 D.x20023_2_0];
% COG = [D.eid 0.774*Reasoning + -0.491*log(ReactionTime)];
COG = [D.eid 0.768*Reasoning + -0.768*log(ReactionTime)];
COG0 = COG(:,1);
COG2 = COG(:,2);

[~, ia] = intersect(COG(:,1), anhedonia_eid(:,1)); cogs{1} = COG(ia, 2); %used time point 0 because it had more data
[~, im] = intersect(COG(:,1), mood_eid(:,1)); cogs{2} = COG(im, 2);
[~, is] = intersect(COG(:,1), somatic_eid(:,1)); cogs{3} = COG(is, 2);
[~, ic] = intersect(COG(:,1), chronic_eid(:,1)); cogs{4} = COG(ic, 2);
[~, il] = intersect(COG(:,1), lateonset_eid(:,1)); cogs{5} = COG(il, 2);
[~, ise] = intersect(COG(:,1), severe_eid(:,1)); cogs{6} = COG(ise, 2);
[~, ih, jh] = intersect(COG(:,1), hetero_eid(:,1)); cogs{7} = COG(ih, 2);
for n = 1:6
    nans{n} = isnan(cogs{n});
end
for n = 1:6
    cog_nanless{n} = cogs{n}(~any(nans{n},2),:);
    kms_cog{n} = kms_final{n}(~any(nans{n},2),:);
end
%run the anova
for i=1:6
        [p_cogs_6(i), C_tbl{i}, stats] = anova1(cog_nanless{i}, kms_cog{i},  'display', 'off'); %this is exactly a ttest, I just did it this way cause I didn't want to restructure my data. To calculate t stat from F stat, t2 = F
F(i) = C_tbl{i}(2,6);
c{i}(:,:) = multcompare(stats);
end
[pthr,pcor,padj] = fdr(p_cogs_6); %%%% fdr.m by Anderson M. Winkler
% Research Imaging Center/UTHSCSA
% Dec/2007 (first version)
% Nov/2012 (this version)
% http://brainder.org

%%%%groups 4 and 6 are significantly different

%% neuroticism (anhedonia was significant but BARELY, prolly not past correction)
D = readtable('kayla.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
N12 = sum([D.x1920_0_0, D.x1930_0_0, D.x1940_0_0, D.x1950_0_0, D.x1960_0_0, D.x1970_0_0, D.x1980_0_0, D.x1990_0_0, D.x2000_0_0, D.x2010_0_0, D.x2020_0_0, D.x2030_0_0], 2);

[~, ia] = intersect(D.eid, anhedonia_eid(:,1)); N{1} = N12(ia, :);
[~, im] = intersect(D.eid, mood_eid(:,1)); N{2} = N12(im, :);
[~, is] = intersect(D.eid, somatic_eid(:,1)); N{3} = N12(is, :);
[~, ic] = intersect(D.eid, chronic_eid(:,1)); N{4} = N12(ic, :);
[~, il] = intersect(D.eid, lateonset_eid(:,1)); N{5} = N12(il, :);
[~, ise] = intersect(D.eid, severe_eid(:,1)); N{6} = N12(ise, :);
%remove missing data
for n = 1:6
    nans{n} = isnan(N{n});
end
for n = 1:6
    N_nanless{n} = N{n}(~any(nans{n},2),:);
    N_zscore{n} = zscore(N_nanless{n});
    kms_N{n} = kms_final{n}(~any(nans{n},2),:);
end
%run anova
for i=1:6
        [p_N_6(i), N_tbl{i}] = anova1(N_zscore{i}, kms_N{i},  'display', 'off');
end

%% townsend
%load kms_final ('finalPCAscores_kms_allage_75jul8.mat')
D = readtable('kayla.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);

[~, ia] = intersect(D.eid, anhedonia_eid(:,1)); T{1} = D.x189_0_0(ia, :);
[~, im] = intersect(D.eid, mood_eid(:,1)); T{2} = D.x189_0_0(im, :);
[~, is] = intersect(D.eid, somatic_eid(:,1)); T{3} = D.x189_0_0(is, :);
[~, ic] = intersect(D.eid, chronic_eid(:,1)); T{4} = D.x189_0_0(ic, :);
[~, il] = intersect(D.eid, lateonset_eid(:,1)); T{5} = D.x189_0_0(il, :);
[~, ise] = intersect(D.eid, severe_eid(:,1)); T{6} = D.x189_0_0(ise, :);

%remove missing data
for n = 1:6
    nans{n} = isnan(T{n});
end
for n = 1:6
    T_nanless{n} = T{n}(~any(nans{n},2),:);
    T_zscore{n} = zscore(T_nanless{n});
    kms_T{n} = kms_final{n}(~any(nans{n},2),:);
end
%run anova
for i=1:6
        [p_T_6(i), T_tbl{i}] = anova1(T_zscore{i}, kms_T{i},  'display', 'off');
end

%% phenotype Ts from Fs
for i = 1:6
   C_t(i) = sqrt(cell2mat(C_tbl{i}(2,5)));
   N_t(i) = sqrt(cell2mat(N_tbl{i}(2,5)));
   T_t(i) = sqrt(cell2mat(T_tbl{i}(2,5)));
end

%% Functions
function [ARI_dist] = pairwise_ARI(list_of_iters,data_per_iter)

n_iters = length(list_of_iters);

ARI_mtx = zeros(n_iters);

for n = 1:7
for i = 1:n_iters
	idx_1 = list_of_iters{i,n};
	data_1 = data_per_iter{i,n};

	for j = (i+1):n_iters
		idx_2 = list_of_iters{j,n};
		data_2 = data_per_iter{i,n};
		[idx_3, ia,ib] = intersect(idx_1,idx_2);

		ARI_val = rand_index(data_1(ia),data_2(ib), 'adjusted');
		ARI_mtx(i,j) = ARI_val;
	end
end

ARI_dist{n} = triu(ARI_mtx,1);
end
end

function [idx_kms, a] = clusk(Z, varargin)
for n = 1:varargin{1}
    kms(:,n) = kmeans(Z, n);
    kms_sil(:,n) = mean(silhouette(Z, kms(:,n)));
end
[~, a] = max(kms_sil);
idx_kms = kms(:,a);
end