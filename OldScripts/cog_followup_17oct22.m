data = 'Cognition.tsv';
phenotype = Reasoning(:,1);
reasoning_p_vals = pheno_valid(data, phenotype);
[~,pcor_reason,~] = fdr(reasoning_p_vals);

%%
phenotype = log(ReactionTime(:,1));
reaction_p_vals = pheno_valid(data, phenotype);
[~,pcor_react,~] = fdr(reaction_p_vals);


function p_P_6 = pheno_valid(data, phenotype)
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data/allage')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Code')
anhedonia_eid = table2array(readtable('HomogeneousGroups_anhedonia.csv'));
mood_eid = table2array(readtable('HomogeneousGroups_mood.csv'));
somatic_eid = table2array(readtable('HomogeneousGroups_somatic.csv'));
chronic_eid = table2array(readtable('HomogeneousGroups_chronic.csv'));
lateonset_eid = table2array(readtable('HomogeneousGroups_lateonset.csv'));
severe_eid = table2array(readtable('HomogeneousGroups_severe.csv'));
hetero_eid = table2array(readtable('HomogeneousGroups_HeterogeneousComparison.csv'));
load('kms_final_17oct22.mat')
D = readtable(data,'FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);

[~, ia] = intersect(D.eid, anhedonia_eid(:,1)); P{1} = phenotype(ia, :);
[~, im] = intersect(D.eid, mood_eid(:,1)); P{2} = phenotype(im, :);
[~, is] = intersect(D.eid, somatic_eid(:,1)); P{3} = phenotype(is, :);
[~, ic] = intersect(D.eid, chronic_eid(:,1)); P{4} = phenotype(ic, :);
[~, il] = intersect(D.eid, lateonset_eid(:,1)); P{5} = phenotype(il, :);
[~, ise] = intersect(D.eid, severe_eid(:,1)); P{6} = phenotype(ise, :);
%remove missing data
for n = 1:6
    nans{n} = isnan(P{n});
end
for n = 1:6
    P_nanless{n} = P{n}(~any(nans{n},2),:);
    P_zscore{n} = zscore(P_nanless{n});
    kms_P{n} = kms_final{n}(~any(nans{n},2),:);
end
%run anova
for i=1:6
        p_P_6(i) = anova1(P_zscore{i}, kms_P{i},  'display', 'off');
end
end