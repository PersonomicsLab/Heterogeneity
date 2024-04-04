%% select only the subjects with imaging data who meet inclusion criteria
load('imagingdata.mat', 'Zs' , 'Zeids') %those with imaging data
%normative deviations are referred to as Z here because they are z-scored

load('include.mat') %all UKB subjects who met the inclusion criteria

%find the intersection
new_eids = intersect(include, subs);
[~, include_idx] = intersect(Zeids, new_eids);
Zs = Zs(include_idx, :);
Zeids = Zeids(include_idx, :);
%% load the ids for each group
anhedonia_eid = table2array(readtable('anhedonia.csv'));
lowmood_eid = table2array(readtable('lowmood.csv'));
somatic_eid = table2array(readtable('somatic.csv'));
chronic_eid = table2array(readtable('chronic.csv'));
lateonset_eid = table2array(readtable('lateonset.csv'));
acuteimpair_eid = table2array(readtable('acuteimpairment.csv'));
heterogeneous_eid = table2array(readtable('heterogeneouscomparison.csv'));

%get the intersection of the complete eid list and eid list for the
%specific group. This also creates a cell array of the subgroup ids
[groups{1}, ia] = intersect(Zeids(:,1), anhedonia_eid(:,1)); Z{1} = Zs(ia, :);
[groups{2}, im] = intersect(Zeids(:,1), lowmood_eid(:,1)); Z{2} = Zs(im, :);
[groups{3}, is] = intersect(Zeids(:,1), somatic_eid(:,1)); Z{3} = Zs(is, :);
[groups{4}, ic] = intersect(Zeids(:,1), chronic_eid(:,1)); Z{4} = Zs(ic, :);
[groups{5}, il] = intersect(Zeids(:,1), lateonset_eid(:,1)); Z{5} = Zs(il, :);
[groups{6}, iai] = intersect(Zeids(:,1), acuteimpair_eid(:,1)); Z{6} = Zs(iai, :);
[~, ih] = intersect(Zeids(:,1), heterogeneous_eid(:,1)); Z{7} = Zs(ih, :);
%% aim 1 anova
% create an array that provides the group number for each subject
G(ia, 1) = 1;
G(im, 1) = 2;
G(is, 1) = 3;
G(ic, 1) = 4;
G(il, 1) = 5;
G(iai, 1) = 6;
G(ih, 1) = 7;
for n = 1:size(Zs,2)
    [anovaPs(n), tbl{n}, stats{n}] = anovan(Zs(:,n),G,'display','off');
    fstat_true(n) = cell2mat(tbl{n}(2,6)); %this will be used for a supplemental analysis
end
save('fstat_true.mat', 'fstat_true')
[p_fdr] = mafdr(anovaPs); %%False Discovery Rate correction
sigidx = find(p_fdr < 0.05); %% index for significant Zs

% determine the significant post hoc comparisons
for i = 1:length(sigidx)
 c{i} = multcompare(stats{sigidx(i)});
 h{i}= find(c{i}(:,6) < 0.0500);
end
for i = 1:length(sigidx)
   sig_groups_posthoc{i} = c{i}(h{i}, 1:2);
   sig_ps{i} = c{i}(h{i},6);
end
%% aim 2 data driven clustering
for i = 1:1000 %number of bootstraps
    for n = 1:6 %number of clinically dissociated groups
        ix = randi(length(Z{n}), [length(Z{n}),1]);
        data = Z{n}(ix,:);
        [~, score, ~, ~, explained]= pca(Z{n});
        esum = zeros(size(Z{n},2),1);
        esum(1) = explained(1);
        for b = 2:size(Z{n},2)
            esum(b) = esum(b-1) + explained(b);
        end
        f = find(esum > 75); %select the number of PCs that explain 75% variance
        fs(i) = f(1,1);
        [~, ka(i,n)] = clusk(score(:,1:f(1,1)), 10); %kmeans clustering, output is best k for the bootstrap
        i
    end
end
%%find optimal cluster number
for i = 1:6
    cluster(i) = mode(ka(:,i));
end

% final kmeans solutions
for n = 1:6
     [~, score_final{n}, ~, ~, explained]= pca(Z{n});
        esum = zeros(size(Z{n},2),1);
        esum(1) = explained(1);
        for b = 2:size(Z{n},2)
            esum(b) = esum(b-1) + explained(b);
        end
        f = find(esum > 75);
        fs(i) = f(1,1);
        [kms_final{n}] = kmeans(score_final{n}(:,1:f(1,1)), cluster(n)); 
end
%% cross validation to determine cluster stability
% split into test and training
folds = 10; %number of cross validation folds
for x = 1:100 %repeated 100 times
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
            esum = zeros(size(Z{n},2),1);
            esum(1) = explained(1);
            for b = 2:size(Z{n},2)
                esum(b) = esum(b-1) + explained(b);
            end
            f = find(esum > 75);
            [kms_temp{x,i,n}, centroid{x,i,n}] = kmeans(score_temp{x,i,n}(:,1:f(1,1)), cluster(n)); 
            score = Z_New{x,i, n}*coeff{x,i,n}; %apply PCA loadings to heldout set
            score_new{x,i,n} = score(:,1:f(1,1));
        end
    end
    
    for i = 1:folds
        for n=1:6
            for j = 1:size(centroid{x,i,n}, 1)
                for k = 1:size(score_new{x,i,n},1)
                    A{x,i,n}{j,k} = sum((centroid{x,i,n}(j,:)-score_new{x,i,n}(k,:)).^2,2); %determine cluster assignment for heldout data
                    [~,cl_idx{x,i,n}(k)] = min(cell2mat(A{x,i,n}(:,k)));
                end
            end
            total{x,i, n} = zeros(1,length(idxTrain{x,i,n}));
            total{x,i, n}(:,idxTrain{x,i,n}) = kms_temp{x,i,n};
            total{x,i, n}(:,idxNew{x,i,n}) = cl_idx{x,i,n};     
        end
    end
end

%%%calculate ARI (adjusted rand index) for folds
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
%% cluster validation of phenotype: cognitive ability
D = readtable('Cognition.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
Reasoning = [D.x20016_0_0];
ReactionTime = [D.x20023_0_0];
COG = [D.eid 0.768*Reasoning + -0.768*log(ReactionTime)]; %as determined in Lyall et al 2016
clear D Reasoning ReactionTime

%%%%%%%%%% select the cog data for each clinically dissociated group
[~, ia] = intersect(COG(:,1), groups{1}(:,1)); cogs{1} = COG(ia, 2); 
[~, im] = intersect(COG(:,1), groups{2}(:,1)); cogs{2} = COG(im, 2);
[~, is] = intersect(COG(:,1), groups{3}(:,1)); cogs{3} = COG(is, 2);
[~, ic] = intersect(COG(:,1), groups{4}(:,1)); cogs{4} = COG(ic, 2);
[~, il] = intersect(COG(:,1), groups{5}(:,1)); cogs{5} = COG(il, 2);
[~, iai] = intersect(COG(:,1), groups{6}(:,1)); cogs{6} = COG(iai, 2);

%remove the missing data
for n = 1:6
    nans{n} = isnan(cogs{n}); 
end
for n = 1:6
    cog_nanless{n} = cogs{n}(~any(nans{n},2),:);
    kms_cog{n} = kms_final{n}(~any(nans{n},2),:); %remove the missing subjects in cluster assignment
end

%%run ttest
for i=1:6
    [~, p_cogs(i), ci, stats] = ttest2(cog_nanless{i}(kms_cog{i}==1), cog_nanless{i}(kms_cog{i}==2));
    tstat_cog(i) = stats.tstat;
end
pfdr_cog = mafdr(p_cogs);
%% cog posthoc - reasoning or reactiontime
D = readtable('Cognition.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
Reasoning = [D.x20016_0_0];
ReactionTime = [D.x20023_0_0];

%%%%%%%%%% select the Reasoning data for each clinically dissociated group
[~, ia] = intersect(D.eid, groups{1}); Reason{1} = Reasoning(ia, 1);
[~, im] = intersect(D.eid, groups{2}); Reason{2} = Reasoning(im, 1);
[~, is] = intersect(D.eid, groups{3}); Reason{3} = Reasoning(is, 1);
[~, ic] = intersect(D.eid, groups{4}); Reason{4} = Reasoning(ic, 1);
[~, il] = intersect(D.eid, groups{5}); Reason{5} = Reasoning(il, 1);
[~, iai] = intersect(D.eid, groups{6}); Reason{6} = Reasoning(iai, 1);

%remove nans
for n = 1:6
    nans{n} = isnan(Reason{n});
end
for n = 1:6
    Reason_nanless{n} = Reason{n}(~any(nans{n},2),:);
    kms_Reason{n} = kms_final{n}(~any(nans{n},2),:);
end

%%run ttest
for i=1:6
    [~, p_Reason_6(i), ~, stat] = ttest2(Reason_nanless{i}(kms_Reason{i}==1), Reason_nanless{i}(kms_Reason{i}==2));
    tstat_Reason(i) = stat.tstat;
end
clear Reason Reason_nanless kms_Reason Reasoning

%%%%%%%%%% select the reactiontime data for each clinically dissociated group
[~, ia] = intersect(D.eid, groups{1}); React{1} = ReactionTime(ia, 1);
[~, im] = intersect(D.eid, groups{2}); React{2} = ReactionTime(im, 1);
[~, is] = intersect(D.eid, groups{3}); React{3} = ReactionTime(is, 1);
[~, ic] = intersect(D.eid, groups{4}); React{4} = ReactionTime(ic, 1);
[~, il] = intersect(D.eid, groups{5}); React{5} = ReactionTime(il, 1);
[~, iai] = intersect(D.eid, groups{6}); React{6} = ReactionTime(iai, 1);
clear D

%remova nans
for n = 1:6
    nans{n} = isnan(React{n});
end
for n = 1:6
    React_nanless{n} = React{n}(~any(nans{n},2),:);
    kms_React{n} = kms_final{n}(~any(nans{n},2),:);
end

%%run ttest
for i=1:6
        [~, p_React_6(i)] = ttest2(React_nanless{i}(kms_React{i}==1), React_nanless{i}(kms_React{i}==2)); 
end
clear ReactionTime React React_nanless kms_React

pfdr_Reason = fdr(p_Reason_6);
pfdr_React = fdr(p_React_6);
%% cluster validation of phenotype: neuroticism
D = readtable('neuroticism.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
N12 = sum([D.x1920_0_0, D.x1930_0_0, D.x1940_0_0, D.x1950_0_0, D.x1960_0_0, D.x1970_0_0, D.x1980_0_0, D.x1990_0_0, D.x2000_0_0, D.x2010_0_0, D.x2020_0_0, D.x2030_0_0], 2);

%%%%%%%%%% select the neuroticism data for each clinically dissociated group
[~, ia] = intersect(D.eid, groups{1}(:,1)); N{1} = N12(ia, :);
[~, im] = intersect(D.eid, groups{2}(:,1)); N{2} = N12(im, :);
[~, is] = intersect(D.eid, groups{3}(:,1)); N{3} = N12(is, :);
[~, ic] = intersect(D.eid, groups{4}(:,1)); N{4} = N12(ic, :);
[~, il] = intersect(D.eid, groups{5}(:,1)); N{5} = N12(il, :);
[~, iai] = intersect(D.eid, groups{6}(:,1)); N{6} = N12(iai, :);

%remove missing data
for n = 1:6
    nans{n} = isnan(N{n});
end
for n = 1:6
    N_nanless{n} = N{n}(~any(nans{n},2),:);
    kms_N{n} = kms_final{n}(~any(nans{n},2),:);
end

%run ttest
for i=1:6
    [~, pR_N_6(i), ~, stats] = ttest2(N_nanless{i}(kms_N{i}==1), N_nanless{i}(kms_N{i}==2)); 
    tstat_N(i) = stats.tstat;
end
pfdr_N = mafdr(p_N_6); 
clear D N12 N N_nanless N_zscore kms_N
%% cluster validation of phenotype: Townsend socioeconomic deprivation
D = readtable('townsend.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);

%%%%%%%%%% select the townsend data for each clinically dissociated group
[~, ia] = intersect(D.eid, groups{1}(:,1)); T{1} = D.x189_0_0(ia, :);
[~, im] = intersect(D.eid, groups{2}(:,1)); T{2} = D.x189_0_0(im, :);
[~, is] = intersect(D.eid, groups{3}(:,1)); T{3} = D.x189_0_0(is, :);
[~, ic] = intersect(D.eid, groups{4}(:,1)); T{4} = D.x189_0_0(ic, :);
[~, il] = intersect(D.eid, groups{5}(:,1)); T{5} = D.x189_0_0(il, :);
[~, iai] = intersect(D.eid, groups{6}(:,1)); T{6} = D.x189_0_0(iai, :);

%remove missing data
for n = 1:6
    nans{n} = isnan(T{n});
end
for n = 1:6
    T_nanless{n} = T{n}(~any(nans{n},2),:);
    kms_T{n} = kms_final{n}(~any(nans{n},2),:);
end
%run ttest
for i=1:6
    [~, p_T_6(i), ~, stats] = ttest2(T_nanless{i}(kms_T{i}==1), T_nanless{i}(kms_T{i}==2)); 
    tstat_T(i) = stats.tstat;
end
pfdr_T = mafdr(p_T_6);
clear D T T_nanless T_zscore kms_T
%% Functions
function [idx_kms, a] = clusk(Z, varargin)
for n = 1:varargin{1}
    kms(:,n) = kmeans(Z, n);
    kms_sil(:,n) = mean(silhouette(Z, kms(:,n)));
end
[~, a] = max(kms_sil);
idx_kms = kms(:,a);
end