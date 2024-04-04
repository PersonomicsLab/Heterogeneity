%% supplement analyses:
%% normative model accuracy measures
Rho = table2array(readtable('Rho_estimate.txt'));
pRho = table2array(readtable('pRho_estimate.txt'));
EXPV= table2array(readtable('EXPV_estimate.txt'));
MSLL = table2array(readtable('MSLL_estimate.txt'));
SMSE = table2array(readtable('SMSE_estimate.txt'));
RMSE = table2array(readtable('RMSE_estimate.txt'));

figure
tiledlayout(3,2)
nexttile
histogram(Rho, 20)
title('Pearson Correlation')
nexttile
histogram(pRho, 20)
title('Pearson Correlation p-value')
nexttile
histogram(EXPV, 20)
title('Explained Variance')
nexttile
histogram(MSLL, 20)
title('Mean Standardized Log Loss')
nexttile
histogram(SMSE, 20)
title('Standardized Mean Squared Error')
nexttile
histogram(RMSE, 20)
title('Root Mean Squared Error')
% print(gcf, 'normativemodeling_accuracy', '-dpng', '-r300')
%% demographics for clinical groups
clear all; close all;

% Load all data
D = readtable('clinicalcovariates.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);

% load eids
load('mainanalysis_results.mat', 'new_eids');
healthy_eids = table2array(readtable('healthycontrols.csv'));
heterogeneous_eid = table2array(readtable('heterogeneouscomparison.csv')); 
[~, iD] = intersect(D.eid, new_eids);
[~, iH, iHealthy] = intersect(D.eid, healthy_eids); 
Dkeep = D(iD, :);
[~, iC] = setdiff(Dkeep.eid, heterogeneous_eid);
Dhealthy = D(iH, :);
Dclinical = Dkeep(iC, :);
clear D

% Measures of interest
RDS = sum([Dkeep.x2050_2_0 Dkeep.x2060_2_0 Dkeep.x2070_2_0 Dkeep.x2080_2_0],2);
NumEp = zeros(size(Dkeep,1),2);
NumEp(:,1) = Dkeep.x4620_2_0; % depression episodes
NumEp(:,2) = Dkeep.x5386_2_0; % unenthusiastic/disinterested episodes
NumEp = max(NumEp,[],2);
NumEp(NumEp>100) = 100;
NumEp(isnan(NumEp)) = 0;
RDSmood = Dkeep.x2050_2_0;
RDSanhedonia = Dkeep.x2060_2_0;
RDSrestless = Dkeep.x2070_2_0;
RDSlethargy = Dkeep.x2080_2_0;

% Confounds
sex = Dkeep.x31_0_0;
age = Dkeep.x21003_2_0;
sex_H = Dhealthy.x31_0_0;
age_H = Dhealthy.x21003_2_0;

sex_C = Dclinical.x31_0_0;
age_C = Dclinical.x21003_2_0;


anhedonia_eid = table2array(readtable('anhedonia.csv')); [~, Ianhedonia] = intersect(Dkeep.eid, anhedonia_eid); 
mood_eid = table2array(readtable('lowmood.csv')); [~, Imood] = intersect(Dkeep.eid, mood_eid); 
somatic_eid = table2array(readtable('somatic.csv')); [~, Isomatic] = intersect(Dkeep.eid, somatic_eid); 
chronic_eid = table2array(readtable('chronic.csv')); [~, Ichronic] = intersect(Dkeep.eid, chronic_eid); 
lateonset_eid = table2array(readtable('lateonset.csv')); [~, Ilateonset] = intersect(Dkeep.eid, lateonset_eid); 
acuteimpair_eid = table2array(readtable('acuteimpairment.csv')); [~, Iacuteimpair] = intersect(Dkeep.eid, acuteimpair_eid); 
heterogeneous_eid = table2array(readtable('heterogeneouscomparison.csv')); [~, Iheterogeneous] = intersect(Dkeep.eid, heterogeneous_eid); 

% Create summary table
summary = table([length(Ianhedonia); length(Imood); length(Isomatic); length(Ichronic); length(Ilateonset); length(Iacuteimpair); length(Iheterogeneous); length(iH)],'RowNames',{'Anhedonia','Low Mood','Somatic','Chronic','Late Onset','Acute Impairment','Heterogeneous comparison', 'Healthy Control'},'VariableNames',{'Number_of_subjects'});  
summary.Mean_age = [mean(age(Ianhedonia)); mean(age(Imood)); mean(age(Isomatic)); mean(age(Ichronic)); mean(age(Ilateonset)); mean(age(Iacuteimpair)); mean(age(Iheterogeneous)); mean(age_H(iHealthy))];
summary.std_age = [std(age(Ianhedonia)); std(age(Imood)); std(age(Isomatic)); std(age(Ichronic)); std(age(Ilateonset)); std(age(Iacuteimpair)); std(age(Iheterogeneous)); std(age_H(iHealthy))];
summary.Sex_male_percent = [sum(sex(Ianhedonia)); sum(sex(Imood)); sum(sex(Isomatic)); sum(sex(Ichronic)); sum(sex(Ilateonset)); sum(sex(Ianhedonia)); sum(sex(Iheterogeneous)); sum(sex_H(iHealthy))];
summary.Sex_male_percent = summary.Sex_male_percent./summary.Number_of_subjects.*100;
summary.Mean_RDSsum = [mean(RDS(Ianhedonia)); mean(RDS(Imood)); mean(RDS(Isomatic)); mean(RDS(Ichronic)); mean(RDS(Ilateonset)); mean(RDS(Iacuteimpair)); mean(RDS(Iheterogeneous)); 4];
summary.Mean_RDSanhedonia = [mean(RDSanhedonia(Ianhedonia)); mean(RDSanhedonia(Imood)); mean(RDSanhedonia(Isomatic)); mean(RDSanhedonia(Ichronic)); mean(RDSanhedonia(Ilateonset)); mean(RDSanhedonia(Iacuteimpair)); mean(RDSanhedonia(Iheterogeneous)); 1];
summary.Mean_RDSmood = [mean(RDSmood(Ianhedonia)); mean(RDSmood(Imood)); mean(RDSmood(Isomatic)); mean(RDSmood(Ichronic)); mean(RDSmood(Ilateonset)); mean(RDSmood(Iacuteimpair)); mean(RDSmood(Iheterogeneous)); 1];
summary.Mean_RDSrestless = [mean(RDSrestless(Ianhedonia)); mean(RDSrestless(Imood)); mean(RDSrestless(Isomatic)); mean(RDSrestless(Ichronic)); mean(RDSrestless(Ilateonset)); mean(RDSrestless(Iacuteimpair)); mean(RDSrestless(Iheterogeneous)); 1];
summary.Mean_RDSlethargy = [mean(RDSlethargy(Ianhedonia)); mean(RDSlethargy(Imood)); mean(RDSlethargy(Isomatic)); mean(RDSlethargy(Ichronic)); mean(RDSlethargy(Ilateonset)); mean(RDSlethargy(Iacuteimpair)); mean(RDSlethargy(Iheterogeneous)); 1];
summary.Mean_AgeOnset = [mean(AgeOnset(Ianhedonia), "omitnan"); mean(AgeOnset(Imood), "omitnan"); mean(AgeOnset(Isomatic), "omitnan"); mean(AgeOnset(Ichronic), "omitnan"); mean(AgeOnset(Ilateonset), "omitnan"); mean(AgeOnset(Iacuteimpair), "omitnan"); mean(AgeOnset(Iheterogeneous), "omitnan"); NaN];
summary.Mean_NumEp = [mean(NumEp(Ianhedonia)); mean(NumEp(Imood)); mean(NumEp(Isomatic)); mean(NumEp(Ichronic)); mean(NumEp(Ilateonset)); mean(NumEp(Iacuteimpair)); mean(NumEp(Iheterogeneous)); 0];
%% Clinical characteristics of dissociated groups and post hoc comparison table
% run mental_health_descriptors_JAMA.m
run('mental_health_descriptors_JAMA')
%% Longitudinal stability of clinical groups 
clear all; close all;

%load necessary data
D = readtable('clinicalcovariates.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
load('mainanalysis_results.mat', 'groups');

%get RDS for all instances
for i = 1:length(groups)
    % extract the RDS scores for each subject in the clinical groups for each instance 
    [~, ia] = intersect(D.eid, groups{i});
    %extract RDS instance 0, 2, and 3
    RDS_i023{i} = D(ia, ismember(D.Properties.VariableNames, {'x2050_0_0', 'x2050_2_0', 'x2050_3_0', 'x2060_0_0', 'x2060_2_0', 'x2060_3_0', 'x2070_0_0', 'x2070_2_0', 'x2070_3_0', 'x2080_0_0', 'x2080_2_0', 'x2080_3_0'}));
end
clear D

%%%%%%% create grouping variables for anova (clinical group and instance)
group = [{'Anhedonia'}; {'Depressed Mood'}; {'Somatic Disturbance'}; {'Chronic'}; {'Late Onset'}; {'Acute Impairment'}];
count = 1;
for j = 1:6
    for i = 1:size(table2array(RDS_i023{j}),1)
        clinicalgroup{count,1} = group{j};
        count = count + 1;
    end
end
clinicalgroup = [ones(size(RDS_i023{1}.x2060_0_0)); 2*ones(size(RDS_i023{2}.x2060_0_0)); 3*ones(size(RDS_i023{3}.x2060_0_0)); 4*ones(size(RDS_i023{4}.x2060_0_0)); 5*ones(size(RDS_i023{5}.x2060_0_0)); 6*ones(size(RDS_i023{6}.x2060_0_0))]; 
clinicalgroupall = [clinicalgroup; clinicalgroup];
instanceall = [ones(size(clinicalgroup)); 2*ones(size(clinicalgroup))];

%select RDS anhedonia
RDSanhedonia = [RDS_i023{1}.x2060_0_0; RDS_i023{2}.x2060_0_0; RDS_i023{3}.x2060_0_0; RDS_i023{4}.x2060_0_0; RDS_i023{5}.x2060_0_0; RDS_i023{6}.x2060_0_0; RDS_i023{1}.x2060_3_0; RDS_i023{2}.x2060_3_0; RDS_i023{3}.x2060_3_0; RDS_i023{4}.x2060_3_0; RDS_i023{5}.x2060_3_0; RDS_i023{6}.x2060_3_0];
%anova for RDS anhedonia
[~, ~, stats] = anovan(RDSanhedonia, {clinicalgroupall, instanceall}, 'model',2, 'varnames',{'Clinical Group','Instance'});
c_anhedonia = multcompare(stats);
% print(gcf, 'anhedonia_longitudinal_stability', '-dpng', '-r300') 
%the group labels will need to be specified by editing the axes properties

RDSlowmood = [RDS_i023{1}.x2050_0_0; RDS_i023{2}.x2050_0_0; RDS_i023{3}.x2050_0_0; RDS_i023{4}.x2050_0_0; RDS_i023{5}.x2050_0_0; RDS_i023{6}.x2050_0_0; RDS_i023{1}.x2050_3_0; RDS_i023{2}.x2050_3_0; RDS_i023{3}.x2050_3_0; RDS_i023{4}.x2050_3_0; RDS_i023{5}.x2050_3_0; RDS_i023{6}.x2050_3_0];
clinicalgroup = [ones(size(RDS_i023{1}.x2050_0_0)); 2*ones(size(RDS_i023{2}.x2050_0_0)); 3*ones(size(RDS_i023{3}.x2050_0_0)); 4*ones(size(RDS_i023{4}.x2050_0_0)); 5*ones(size(RDS_i023{5}.x2050_0_0)); 6*ones(size(RDS_i023{6}.x2050_0_0))]; 
clinicalgroupall = [clinicalgroup; clinicalgroup];
instanceall = [ones(size(clinicalgroup)); 2*ones(size(clinicalgroup))];
%anova for RDS depressed mood
[~, ~, stats] = anovan(RDSlowmood, {clinicalgroupall, instanceall}, 'model',2, 'varnames',{'clinical group','instance'});
c_mood = multcompare(stats);
% print(gcf, 'lowmood_longitudinal_stability', '-dpng', '-r300')

%select RDS somatic
RDSsomatic = [RDS_i023{1}.x2070_0_0; RDS_i023{2}.x2070_0_0; RDS_i023{3}.x2070_0_0; RDS_i023{4}.x2070_0_0; RDS_i023{5}.x2070_0_0; RDS_i023{6}.x2070_0_0; RDS_i023{1}.x2070_3_0; RDS_i023{2}.x2070_3_0; RDS_i023{3}.x2070_3_0; RDS_i023{4}.x2070_3_0; RDS_i023{5}.x2070_3_0; RDS_i023{6}.x2070_3_0];
clinicalgroup = [ones(size(RDS_i023{1}.x2070_0_0)); 2*ones(size(RDS_i023{2}.x2070_0_0)); 3*ones(size(RDS_i023{3}.x2070_0_0)); 4*ones(size(RDS_i023{4}.x2070_0_0)); 5*ones(size(RDS_i023{5}.x2070_0_0)); 6*ones(size(RDS_i023{6}.x2070_0_0))]; 
clinicalgroupall = [clinicalgroup; clinicalgroup];
instanceall = [ones(size(clinicalgroup)); 2*ones(size(clinicalgroup))];
%anova for RDS somatic
[~, ~, stats] = anovan(RDSsomatic, {clinicalgroupall, instanceall}, 'model',2, 'random',[1,2], 'varnames',{'clinical group','instance'});
c_somatic = multcompare(stats);
% print(gcf, 'somatic_longitudinal_stability', '-dpng', '-r300')
%% Postdoc anova comparisons of main aim 1 analysis
load('aim1results_21mar24.mat', 'sig_groups_posthoc', 'sig_ps')
%% T-tests against 0 for each clinically dissociated group
% run ttest_against0.m
run('ttest_against0')
%% Reduce heterogeneous group for aim 1 result robustness 
clear all; close all;

% load necessary data
load('imagingdata.mat', 'Z')
load('mainanalysis_results.mat', 'sigidx', 'fignames')

% find average size of groups
for i = 1:6
   count(i) = size(Z{i}, 1);
end
avg = round(sum(count)/6);

% run bootstraps
heterogeneousZs = Z{7};
for i = 1:1000
    for n = 1:size(sigidx,2)
        ix = randi(size(heterogeneousZs, 1), [avg,1]); %300 is the average of 6 groups 
        subsamp = Z{7}(ix, sigidx(n));
        bootmeans(i, n) = mean(subsamp);
    end
        i
end

%plot figure
figure 
tiledlayout(2,5)
for i = 1:length(sigidx)
    nexttile
    histogram(bootmeans(:,i))
    title(fignames{i})
    xlabel('Normative Deviation Mean')
    vline([mean(Z{1}(:,i)) mean(Z{2}(:,i)) mean(Z{3}(:,i)) mean(Z{4}(:,i)) mean(Z{5}(:,i)) mean(Z{6}(:,i))], {'b'}) %,{[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]},{'anhedonia','depressed mood','somatic', 'chronic', 'late onset', 'acute'})
    sorted = sort(bootmeans(:,i), 'descend');
    vline(sorted(25), 'r', 'Top 2.5%')
    vline(sorted(975), 'r', 'Bottom 2.5%')
end
% print(gcf, 'reduce_heterogeneous_N', '-dpng', '-r300')
%% held out demographics
clear all; close all;
D = readtable('clinicalcovariates_forsubgroups_23feb22.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
% load eids
load('heldout_eids.mat')
%keep D for eids
[~, iD] = intersect(D.eid, all_eid_heldout);
Dkeep = D(iD, :);
clear D

% Measures of interest
RDS = sum([Dkeep.x2050_2_0 Dkeep.x2060_2_0 Dkeep.x2070_2_0 Dkeep.x2080_2_0],2);
AgeOnset = Dkeep.x20433_0_0;
NumEp = zeros(size(Dkeep,1),2);
NumEp(:,1) = Dkeep.x4620_2_0; % depression episodes
NumEp(:,2) = Dkeep.x5386_2_0; % unenthusiastic/disinterested episodes
NumEp = max(NumEp,[],2);
NumEp(NumEp>100) = 100;
NumEp(isnan(NumEp)) = 0;
RDSmood = Dkeep.x2050_2_0;
RDSanhedonia = Dkeep.x2060_2_0;
RDSrestless = Dkeep.x2070_2_0;
RDSlethargy = Dkeep.x2080_2_0;

% Confounds
sex = Dkeep.x31_0_0;
age = Dkeep.x21003_2_0;

% find the index for each clinical groupss
[~, Ianhedonia] = intersect(Dkeep.eid, a_heldout_eids); 
[~, Imood] = intersect(Dkeep.eid, m_heldout_eids); 
[~, Isomatic] = intersect(Dkeep.eid, s_heldout_eids); 
[~, Ichronic] = intersect(Dkeep.eid, c_heldout_eids); 
[~, Ilateonset] = intersect(Dkeep.eid, l_heldout_eids); 
[~, Iacuteimpair] = intersect(Dkeep.eid, ai_heldout_eids); 
[~, Iheterogeneous] = intersect(Dkeep.eid, h_heldout_eids); 

% Create summary table
summary = table([length(Ianhedonia); length(Imood); length(Isomatic); length(Ichronic); length(Ilateonset); length(Iacuteimpair); length(Iheterogeneous)],'RowNames',{'Anhedonia','Low Mood','Somatic','Chronic','Late Onset','Acute Impairment','Heterogeneous comparison'},'VariableNames',{'Number_of_subjects'});  
summary.Mean_age = [mean(age(Ianhedonia)); mean(age(Imood)); mean(age(Isomatic)); mean(age(Ichronic)); mean(age(Ilateonset)); mean(age(Iacuteimpair)); mean(age(Iheterogeneous))];
summary.std_age = [std(age(Ianhedonia)); std(age(Imood)); std(age(Isomatic)); std(age(Ichronic)); std(age(Ilateonset)); std(age(Iacuteimpair)); std(age(Iheterogeneous))];
summary.Sex_male_percent = [sum(sex(Ianhedonia)); sum(sex(Imood)); sum(sex(Isomatic)); sum(sex(Ichronic)); sum(sex(Ilateonset)); sum(sex(Ianhedonia)); sum(sex(Iheterogeneous))];
summary.Sex_male_percent = summary.Sex_male_percent./summary.Number_of_subjects.*100;
summary.Mean_RDSsum = [mean(RDS(Ianhedonia)); mean(RDS(Imood)); mean(RDS(Isomatic)); mean(RDS(Ichronic)); mean(RDS(Ilateonset)); mean(RDS(Iacuteimpair)); mean(RDS(Iheterogeneous))];
summary.Mean_RDSanhedonia = [mean(RDSanhedonia(Ianhedonia)); mean(RDSanhedonia(Imood)); mean(RDSanhedonia(Isomatic)); mean(RDSanhedonia(Ichronic)); mean(RDSanhedonia(Ilateonset)); mean(RDSanhedonia(Iacuteimpair)); mean(RDSanhedonia(Iheterogeneous))];
summary.Mean_RDSmood = [mean(RDSmood(Ianhedonia)); mean(RDSmood(Imood)); mean(RDSmood(Isomatic)); mean(RDSmood(Ichronic)); mean(RDSmood(Ilateonset)); mean(RDSmood(Iacuteimpair)); mean(RDSmood(Iheterogeneous))];
summary.Mean_RDSrestless = [mean(RDSrestless(Ianhedonia)); mean(RDSrestless(Imood)); mean(RDSrestless(Isomatic)); mean(RDSrestless(Ichronic)); mean(RDSrestless(Ilateonset)); mean(RDSrestless(Iacuteimpair)); mean(RDSrestless(Iheterogeneous))];
summary.Mean_RDSlethargy = [mean(RDSlethargy(Ianhedonia)); mean(RDSlethargy(Imood)); mean(RDSlethargy(Isomatic)); mean(RDSlethargy(Ichronic)); mean(RDSlethargy(Ilateonset)); mean(RDSlethargy(Iacuteimpair)); mean(RDSlethargy(Iheterogeneous))];
summary.Mean_AgeOnset = [mean(AgeOnset(Ianhedonia), "omitnan"); mean(AgeOnset(Imood), "omitnan"); mean(AgeOnset(Isomatic), "omitnan"); mean(AgeOnset(Ichronic), "omitnan"); mean(AgeOnset(Ilateonset), "omitnan"); mean(AgeOnset(Iacuteimpair), "omitnan"); mean(AgeOnset(Iheterogeneous), "omitnan")];
summary.Mean_NumEp = [mean(NumEp(Ianhedonia)); mean(NumEp(Imood)); mean(NumEp(Isomatic)); mean(NumEp(Ichronic)); mean(NumEp(Ilateonset)); mean(NumEp(Iacuteimpair)); mean(NumEp(Iheterogeneous))];
%% Held out reduce heterogenous group
clear Z Zever subsamp bootmeans count

% load Zs
load('heldout_Z.mat')

%%select Zs for groups
[~, ia] = intersect(all_eid_heldout(:,1), a_heldout_eids(:,1)); Z{1} = Zs(ia, :);
[~, im] = intersect(all_eid_heldout(:,1), m_heldout_eids(:,1)); Z{2} = Zs(im, :);
[~, is] = intersect(all_eid_heldout(:,1), s_heldout_eids(:,1)); Z{3} = Zs(is, :);
[~, ic] = intersect(all_eid_heldout(:,1), c_heldout_eids(:,1)); Z{4} = Zs(ic, :);
[~, il] = intersect(all_eid_heldout(:,1), l_heldout_eids(:,1)); Z{5} = Zs(il, :);
[~, iai] = intersect(all_eid_heldout(:,1), ai_heldout_eids(:,1)); Z{6} = Zs(iai, :);
[~, ih] = intersect(all_eid_heldout(:,1), h_heldout_eids(:,1)); Z{7} = Zs(ih, :);

load('mainanalysis_results.mat', 'sigidx', 'fignames')

for i = 1:1000
    for n = 1:length(sigidx)
        ix = randi(size(heterogeneousZs,1), [size(heterogeneousZs,1), 1]); 
        subsamp = Z{7}(ix, sigidx(n));
        bootmeans(i, n) = mean(subsamp);
    end
        i
end

figure 
tiledlayout(2,5)
for i = 1:length(sigidx)
    nexttile
    histogram(bootmeans(:,i))
    title(fignames{i})
    xlabel('Normative Deviation Mean')
    vline([mean(Z{1}(:,i)) mean(Z{2}(:,i)) mean(Z{3}(:,i)) mean(Z{4}(:,i)) mean(Z{5}(:,i)) mean(Z{6}(:,i))], {'b'})
    sorted = sort(bootmeans(:,i), 'descend');
    vline(sorted(25), 'r', 'Top 2.5%')
    vline(sorted(975), 'r', 'Bottom 2.5%')
end
% print(gcf, 'reduce_heterogeneous_N_heldout', '-dpng', '-r300')
%% ANOVA comparisons against random subgroups
close all; clear all;

%load data
load('imagingdata.mat', 'Zs')
load('mainanalysis_results.mat', 'sigidx', 'fignames')

for i = 1:1000
    ix = randperm(size(Zs,1));
    group = G(ix,:);
    for n = 1:size(sigidx,2)
        [~, tbl] = anovan(Zs(:,sigidx(n)),group,'display','off'); % run ANOVA
        fstat(i,n) = cell2mat(tbl(2,6));
    end
end
%find top 5% f-stat for each of the sig idps
sortedData = sort(fstat, 'descend');
top_threshold = sortedData(25, :);

%load the true fstats
load('fstat_true.mat')

%plot figure
figure
tiledlayout(2,5);
for i = 1:length(sigidx)
    nexttile
    histogram(fstat(:,i))
    title(fignames{i})
    xlabel('F-statistic')
    vline(fstat_true(sigidx(i)))
    vline(top_threshold(i), 'k', '2.5%')
end
%print(gcf,'randomgrouping_fstat_bootstrap','-dpng','-r300');
%% Clustering silhouette scores
clear all; close all;
load('imagingdata.mat', 'Z')
load('mainanalysis_results.mat', 'sigidx')

for i = 1:1000
    for n = 1:6
        ix = randi(size(Z{n},1), [size(Z{n}, 1),1]);
        [~, ka(i,n)] = clusk(Z{n}(ix,sigidx), 10); 
    end
            i
end

%%histogram
figure
hist(ka)
xlabel('Best Cluster Number')
ylabel('Number of Boostraps')
% print(gcf, 'mode_clust_assign_figure', '-dpng', '-r300')
%% cluster spider plots
%%%%%%%% need a finalized cluster assignment
clear all; close all;
load('imagingdata.mat', 'Z')
load('mainanalysis_results.mat', 'kms_final')

for i = 1:6
    meanZ{i}(:,1) = mean(Z{i}(kms_final{i} == 1, :),1);
    meanZ{i}(:,2) = mean(Z{i}(kms_final{i} == 2, :),1);
end

%%seperate by modality
for i = 1:6
    mean_mods{i} = [mean(meanZ{i}(1:20,:),1); mean(meanZ{i}(20:50,:),1); mean(meanZ{i}(51:62,:),1); mean(meanZ{i}(63:90,:),1)];
end
mod_names = [{'FA', 'GMV', 'CT', 'FC'}];

%%visualization, axes between modalities are consistent
group = [{'Anhedonia', 'Low Mood', 'Somatic', 'Chronic', 'Late Onset', 'Acute Impairment'}];
figure
tiledlayout(2,3)
for i = 1:6
    nexttile
    spider_plot(mean_mods{i}(:,1:2)','AxesLabels', mod_names, 'LabelFontSize', 12, 'AxesPrecision', 2, 'AxesLabelsEdge', 'none', 'LabelFontSize', 12, 'AxesLimits', [-.7, -.7, -.7, -.7; .5, .5, .5, .5])
    title(group{i})
end
    legend({'Cluster 1', 'Cluster 2'})
%print(gcf,'aim2_clust_spiderplot', '-dpng','-r300');
%% aim 2 cluster phenotye ttests post hoc comparisons
load('mainanalysis_results.mat', 'tstat_cog', 'pfdr_cog', 'tstat_N', 'pfdr_N', 'tstat_T', 'pfdr_T')
summary = table([tstat_cog; pfdr_cog; tstat_N; pfdr_T; tstat_T; pfdr_T], 'RowNames', {'Anhedonia', 'Low Mood', 'Somatic', 'Chronic', 'Late Onset', 'Acute Impairment'});