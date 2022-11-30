addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Code')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data/allage/')

clear all; close all; clc

% Load all data
%warning off
D = readtable('Subgroups.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
IDP = readtable('IDPs_select89_40k.tsv','FileType','text');
if sum(D.eid - table2array(IDP(:,1))) ~= 0 
    fprintf('WARNING: subject IDs do not line up\n');
end 
subIDs = table2array(IDP(:,1));
IDP = table2array(IDP(:,2:end));
IDP(IDP==0) = nan;

% Measures of interest
RDS = sum([D.x2050_2_0 D.x2060_2_0 D.x2070_2_0 D.x2080_2_0],2);
AgeOnset = D.x20433_0_0;
AgeOnset(isnan(AgeOnset))=0;
NumEp = zeros(size(D,1),2);
NumEp(:,1) = D.x4620_2_0; % depression episodes
NumEp(:,2) = D.x5386_2_0; % unenthusiastic/disinterested episodes
NumEp = max(NumEp,[],2);
NumEp(NumEp>100) = 100;
NumEp(isnan(NumEp)) = 0;
WMH = D.x25781_2_0;
RDSmood = D.x2050_2_0;
RDSanhedonia = D.x2060_2_0;
RDSrestless = D.x2070_2_0;
RDSlethargy = D.x2080_2_0;

% Confounds
sex = D.x31_0_0;
age = D.x21003_2_0;
motion = D.x25741_2_0;
MedOnOff = Match_medication(D.x20003_2_0);

% Remove subjects over 30 and with missing brain data 
Ibrain = find(isnan(sum(IDP,2))==0);
Irds = find(isnan(RDS)==0);
% Ilatelife = find(D.x21003_2_0>=60); 
% Ikeep = intersect(Ibrain, Ilatelife); Ikeep = intersect(Ikeep, Irds);
Iallage = find(D.x21003_2_0>=30);
Ikeep = intersect(Ibrain, Iallage); Ikeep = intersect(Ikeep, Irds);
IDP = IDP(Ikeep,:); subIDs = subIDs(Ikeep,:);
RDS = RDS(Ikeep,:); RDSmood = RDSmood(Ikeep,:); RDSanhedonia = RDSanhedonia(Ikeep,:); RDSlethargy = RDSlethargy(Ikeep,:); RDSrestless = RDSrestless(Ikeep,:);
AgeOnset = AgeOnset(Ikeep,:); NumEp = NumEp(Ikeep,:); WMH = WMH(Ikeep,:); 
MedOnOff = MedOnOff(Ikeep,:); sex = sex(Ikeep,:); age = age(Ikeep,:); motion = motion(Ikeep,:);
clear Ilatelife D Ibrain %Ikeep

% WMH abnormalities
ages = [60 65 70 75 80 100];
WMHextreme = nan(length(WMH),1);
for a = 1:length(ages)-1
    I = find(age>=ages(a) & age<ages(a+1));
    m = nanmean(WMH(I));
    s = nanstd(WMH(I));
    up = find(WMH(I) >= m+2*s);
    down = find(WMH(I) <=m );
    WMHextreme(I(up)) = 1;
    WMHextreme(I(down)) = -1;
    WMHextreme(I(setdiff(1:length(I),[up;down]))) = 0;
end
clear ages a I m s up down

% Extreme subgroups
Ilateonset = find(AgeOnset>=60 & RDS<=5 & NumEp<=10); dlmwrite('HomogeneousGroups_lateonset.csv',subIDs(Ilateonset),'precision',7);
% Iearlyonset = find(AgeOnset<60 & RDS<=5 & NumEp<=10);
Ichronic = find(NumEp>10 & RDS<=5 & AgeOnset<60); dlmwrite('HomogeneousGroups_chronic.csv',subIDs(Ichronic),'precision',7);
Imood = find(RDSmood>2 & RDSrestless<=2 & RDSanhedonia<=2 & AgeOnset<60 & NumEp<=10); dlmwrite('HomogeneousGroups_mood.csv',subIDs(Imood),'precision',7);
Isomatic = find(RDSmood<=2 & RDSrestless>2 & RDSanhedonia<=2 & AgeOnset<60 & NumEp<=10); dlmwrite('HomogeneousGroups_somatic.csv',subIDs(Isomatic),'precision',7);
Ianhedonia = find(RDSmood<=2 & RDSrestless<=2 & RDSanhedonia>2 & AgeOnset<60 & NumEp<=10); dlmwrite('HomogeneousGroups_anhedonia.csv',subIDs(Ianhedonia),'precision',7);
Isevere = find(RDS>=10 & AgeOnset<60 & NumEp<=10); Isevere = setdiff(Isevere,unique([Ilateonset; Ichronic; Imood; Isomatic; Ianhedonia])); dlmwrite('HomogeneousGroups_severe.csv',subIDs(Isevere),'precision',7);
Inormal = find(RDS==4 & AgeOnset==0 & NumEp==0 & MedOnOff==0); dlmwrite('HomogeneousGroups_healthycontrols.csv',subIDs(Inormal),'precision',7);
Iheterogeneous = setdiff(1:length(subIDs),unique([Ilateonset; Ichronic; Imood; Isomatic; Ianhedonia; Isevere; Inormal; find(RDS==4); find(AgeOnset==0); find(NumEp==0)])); dlmwrite('HomogeneousGroups_HeterogeneousComparison.csv',subIDs(Iheterogeneous),'precision',7);
Icombined_for_wei = [Ianhedonia; Imood; Isomatic; Ichronic; Ilateonset; Isevere; Iheterogeneous'];
dlmwrite('HomogeneousGroups_all.csv',[subIDs(Icombined_for_wei) RDS(Icombined_for_wei)],'precision',7);
IDcombined = subIDs(Icombined_for_wei);
IDcombined(:,2) = [ones(length(Ianhedonia),1); 2*ones(length(Imood),1); 3*ones(length(Isomatic),1); 4*ones(length(Ichronic),1); 5*ones(length(Ilateonset),1); 6*ones(length(Isevere),1); 7*ones(length(Iheterogeneous),1)]; 
% anhedonia_eid = table2array(readtable('anhedonia_eid.csv'));
% mood_eid = table2array(readtable('mood_eid.csv'));
% somatic_eid = table2array(readtable('somatic_eid.csv'));
% chronic_eid = table2array(readtable('chronic_eid.csv'));
% lateonset_eid = table2array(readtable('lateonset_eid.csv'));
% severe_eid = table2array(readtable('severe_eid.csv'));
% hetero_eid = table2array(readtable('hetero_eid.csv'));
% IDcombined = [anhedonia_eid(:,2); mood_eid(:,2); somatic_eid(:,2); chronic_eid(:,2); lateonset_eid(:,2); severe_eid(:,2); hetero_eid(:,2)];
% IDcombined(:,2) = [ones(195,1); 2*ones(144,1); 3*ones(179,1); 4*ones(235,1); 5*ones(529,1); 6*ones(307,1); 7*ones(2445,1)]; 
% [IDz,i] = sort(IDcombined(:,1)); G = IDcombined(i,2);

% Find duplicates
ID = subIDs(Icombined_for_wei);
[uniqueA, i, j] = unique(ID,'first');
indexToDupes = find(not(ismember(1:numel(ID),i)))';
Dupes = ID(indexToDupes);
for n = 1:length(Dupes)
    indexToDupes(n,2) = find(ID==Dupes(n),1,'first');
end
Dupes = [Dupes indexToDupes]; clear uniqueA i j n indexToDupes

% Create summary table
summary = table([length(Ilateonset); length(Ichronic); length(Imood); length(Isomatic); length(Ianhedonia); length(Isevere); length(Inormal); length(Iheterogeneous)],'RowNames',{'Late_onset','Chronic','Symptom_mood','Symptom_somatic','Symptom_anhedonia','Acute_severity','HC','Heterogeneous_comparison'},'VariableNames',{'Number_of_subjects'});  
summary.Mean_age = [mean(age(Ilateonset)); mean(age(Ichronic)); mean(age(Imood)); mean(age(Isomatic)); mean(age(Ianhedonia)); mean(age(Isevere)); mean(age(Inormal)); mean(age(Iheterogeneous))];
summary.std_age = [std(age(Ilateonset)); std(age(Ichronic)); std(age(Imood)); std(age(Isomatic)); std(age(Ianhedonia)); std(age(Isevere)); std(age(Inormal)); std(age(Iheterogeneous))];
summary.Sex_male_percent = [sum(sex(Ilateonset)); sum(sex(Ichronic)); sum(sex(Imood)); sum(sex(Isomatic)); sum(sex(Isevere)); sum(sex(Ianhedonia)); sum(sex(Inormal)); sum(sex(Iheterogeneous))];
summary.Sex_male_percent = summary.Sex_male_percent./summary.Number_of_subjects.*100;
summary.Head_motion = [mean(motion(Ilateonset)); mean(motion(Ichronic)); mean(motion(Imood)); mean(motion(Isomatic)); mean(motion(Ianhedonia)); mean(motion(Isevere)); mean(motion(Inormal)); mean(motion(Iheterogeneous))];
summary.WMH = [mean(WMH(Ilateonset)); mean(WMH(Ichronic)); mean(WMH(Imood)); mean(WMH(Isomatic)); mean(WMH(Ianhedonia)); mean(WMH(Isevere)); mean(WMH(Inormal)); mean(WMH(Iheterogeneous))];
summary.Medication_on_percent = [sum(MedOnOff(Ilateonset)); sum(MedOnOff(Ichronic)); sum(MedOnOff(Imood)); sum(MedOnOff(Isomatic)); sum(MedOnOff(Ianhedonia)); sum(MedOnOff(Isevere)); sum(MedOnOff(Inormal)); sum(MedOnOff(Iheterogeneous))];
summary.Medication_on_percent = summary.Medication_on_percent./summary.Number_of_subjects.*100;
summary.Mean_RDS = [mean(RDS(Ilateonset)); mean(RDS(Ichronic)); mean(RDS(Imood)); mean(RDS(Isomatic)); mean(RDS(Ianhedonia)); mean(RDS(Isevere)); mean(RDS(Inormal)); mean(RDS(Iheterogeneous))];


%% summary for clusters
summary = table([length(Ianhedonia(kms_final{1}==1)); length(Ianhedonia(kms_final{1}==2)); length(Imood(kms_final{2}==1)); length(Imood(kms_final{2}==2)); length(Isomatic(kms_final{3}==1)); length(Isomatic(kms_final{3}==2)); length(Ichronic(kms_final{4}==1)); length(Ichronic(kms_final{4}==2)); length(Ilateonset(kms_final{5}==1)); length(Ilateonset(kms_final{5}==2)); length(Isevere(kms_final{6}==1)); length(Isevere(kms_final{6}==2))],'RowNames',{'Symptom_anhedonia C1','Symptom_anhedonia C2','Symptom_mood C1','Symptom_mood C2','Symptom_somatic C1','Symptom_somatic C2','Chronic C1','Chronic C2','Late_Onset C1','Late_Onset C2','Acute_severity C1', 'Acute_severity C2'},'VariableNames',{'Number_of_subjects'});  
summary.Mean_age = [mean(age(Ianhedonia(kms_final{1}==1))); mean(age(Ianhedonia(kms_final{1}==2))); mean(age(Imood(kms_final{2}==1))); mean(age(Imood(kms_final{2}==2))); mean(age(Isomatic(kms_final{3}==1))); mean(age(Isomatic(kms_final{3}==2))); mean(age(Ichronic(kms_final{4}==1))); mean(age(Ichronic(kms_final{4}==2))); mean(age(Ilateonset(kms_final{5}==1))); mean(age(Ilateonset(kms_final{5}==2))); mean(age(Isevere(kms_final{6}==1))); mean(age(Isevere(kms_final{6}==2)))];
summary.Std_age = [std(age(Ianhedonia(kms_final{1}==1))); std(age(Ianhedonia(kms_final{1}==2))); std(age(Imood(kms_final{2}==1))); std(age(Imood(kms_final{2}==2))); std(age(Isomatic(kms_final{3}==1))); std(age(Isomatic(kms_final{3}==2))); std(age(Ichronic(kms_final{4}==1))); std(age(Ichronic(kms_final{4}==2))); std(age(Ilateonset(kms_final{5}==1))); std(age(Ilateonset(kms_final{5}==2))); std(age(Isevere(kms_final{6}==1))); std(age(Isevere(kms_final{6}==1)))];
summary.Sex_male_percent = [sum(sex(Ianhedonia(kms_final{1}==1))); sum(sex(Ianhedonia(kms_final{1}==2))); sum(sex(Imood(kms_final{2}==1))); sum(sex(Imood(kms_final{2}==2))); sum(sex(Isomatic(kms_final{3}==1))); sum(sex(Isomatic(kms_final{3}==2))); sum(sex(Ichronic(kms_final{4}==1))); sum(sex(Ichronic(kms_final{4}==2))); sum(sex(Ilateonset(kms_final{5}==1))); sum(sex(Ilateonset(kms_final{5}==2))); sum(sex(Isevere(kms_final{6}==1))); sum(sex(Isevere(kms_final{6}==1)))];
summary.Sex_male_percent = summary.Sex_male_percent./summary.Number_of_subjects.*100;
summary.Mean_RDS = [mean(RDS(Ianhedonia(kms_final{1}==1))); mean(RDS(Ianhedonia(kms_final{1}==2))); mean(RDS(Imood(kms_final{2}==1))); mean(RDS(Imood(kms_final{2}==2))); mean(RDS(Isomatic(kms_final{3}==1))); mean(RDS(Isomatic(kms_final{3}==2))); mean(RDS(Ichronic(kms_final{4}==1))); mean(RDS(Ichronic(kms_final{4}==2))); mean(RDS(Ilateonset(kms_final{5}==1))); mean(RDS(Ilateonset(kms_final{5}==2))); mean(RDS(Isevere(kms_final{6}==1))); mean(RDS(Isevere(kms_final{6}==2)))];
