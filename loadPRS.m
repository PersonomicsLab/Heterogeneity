% Load PRS
D1 = readtable('/scratch/khannon/TNM_Subject_IDs.csv');
%D1 = subject_eid.eid;
Ryan_MDD = readtable('/scratch/khannon/UKB_PRS/mvp_mdd_ukb_eur_prs_complete.profile','FileType','text');
Ryan_AD = readtable('/scratch/khannon/UKB_PRS/jansen_ad_ukb_eur_prs_complete.profile','FileType','text');
Ryan_fam = readtable('/scratch/khannon/UKB_PRS/ukb19798_cal_chr22_v2_s488265.fam','FileType','text');
Janine_fam = readtable('/ceph/biobank/phenotypes/ukb47267_cal_chr17_v2_s488282.fam','FileType','text');
PRS_MDD = nan(size(D1,1),2); PRS_AD = nan(size(D1,1),2);
for n = 1:size(D1,1)
    I1 = find(Janine_fam.Var1==D1.eid(n));
    if I1
        I2 = find(Ryan_MDD.FID==Ryan_fam.Var1(I1));
        if I2
            PRS_MDD(n,1) = D1.eid(I1);
            PRS_MDD(n,2) = Ryan_MDD.SCORE(I2);
            PRS_AD(n,2) = Ryan_AD.SCORE(I2);
        end
    end
end
I = find(isnan(PRS_AD)==0);
PRS_AD = PRS_AD(I);
PRS_MDD = PRS_MDD(I);
save('May23_PRS.mat')

%%
%%load eid
anhedonia_eid = table2array(readtable('anhedonia_eid.csv'));
mood_eid = table2array(readtable('mood_eid.csv'));
somatic_eid = table2array(readtable('somatic_eid.csv'));
chronic_eid = table2array(readtable('chronic_eid.csv'));
lateonset_eid = table2array(readtable('lateonset_eid.csv'));
severe_eid = table2array(readtable('severe_eid.csv'));
hetero_eid = table2array(readtable('hetero_eid.csv'));

%%sort and find the eids for each group
D1 = readtable('TNM_Subject_IDs.csv');
D = D1.eid(I);
PRS = PRS_MDD;
[~, ia, ja] = intersect(D, anhedonia_eid(:,2)); PRS_anhedonia = PRS(ia, :);
[~, im, jm] = intersect(D, mood_eid(:,2)); PRS_mood = PRS(im, :);
[~, is, js] = intersect(D, somatic_eid(:,2)); PRS_somatic = PRS(is, :);
[~, ic, jc] = intersect(D, chronic_eid(:,2)); PRS_chronic = PRS(ic, :);
[~, il, jl] = intersect(D, lateonset_eid(:,2)); PRS_lateonset = PRS(il, :);
[~, ise, jse] = intersect(D, severe_eid(:,2)); PRS_severe = PRS(ise, :);
% PRS_MDD_anhedonia = PRS_MDD(ia, :);
% PRS_MDD_mood = PRS_MDD(im, :);
% PRS_MDD_somatic = PRS_MDD(is, :);
% PRS_MDD_chronic = PRS_MDD(ic, :);
% PRS_MDD_lateonset = PRS_MDD(il, :);
% PRS_MDD_severe = PRS_MDD(ise, :);

index{1} = ja;
index{2} = jm;
index{3} = js;
index{4} = jc;
index{5} = jl;
index{6} = jse;

PRSMDD{1} = PRS_anhedonia;
PRSMDD{2} = PRS_mood;
PRSMDD{3} = PRS_somatic;
PRSMDD{4} = PRS_chronic;
PRSMDD{5} = PRS_lateonset;
PRSMDD{6} = PRS_severe;

% PRSMDD{1} = PRS_MDD_anhedonia;
% PRSMDD{2} = PRS_MDD_mood;
% PRSMDD{3} = PRS_MDD_somatic;
% PRSMDD{4} = PRS_MDD_chronic;
% PRSMDD{5} = PRS_MDD_lateonset;
% PRSMDD{6} = PRS_MDD_severe;
%%working on how to select only kms of the data we have
%kms is a cell so I need to make the eids cells
%index iterates through the whoooole thing, not just the group uuuughh
%how do I make the index specific to the cell format? I'm not sure that's
%even doable
%maybe it's more worthwhile to move kms to long format instead of the other
%way around - hm I don't think I can. 
%okay so don't apply index onto kms, apply it onto something else first to
%get an index in the right format
%okay so I have an index of the long format, how do I convert it to short
%format?
for i=1:6
    kms_PRSMDD{i} = kms_TNM_3_all{i}(index{i});
end
%shit that's not good, why is index iterating past the 6 groups
%%
km = [kms{1}; kms{2}; kms{3}; kms{4}; kms{5}; kms{6}];
%%
% for n = 1:6
%     kms_nanless{n} = kms{n}(~any(nans{n},2),:);
% end


%%anova
%this anove code is not right lol but I don't know what's wrong, nans?
%okay I need to make sure the anova is properly recognizing that there are
%multiple clsuters to compare within each group
%okay so kms_nanless needs to be a 
for i=1:6
        p_PRSAD_6(i) = anova1(PRS{i}, kms_B{i},  'display', 'off');
end
for i=1:6
        p_PRSMDD_6(i) = anova1(PRSMDD{i}, kms_B{i},  'display', 'off');
end

%% cognition
%load cognition again. select the 2, multiply the 2 by the weights
D = readtable('Cognition.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
Reasoning = [D.x20016_0_0 D.x20016_2_0];
ReactionTime = [D.x20023_0_0 D.x20023_2_0];
COG = [D.eid 0.774*Reasoning + -0.491*log(ReactionTime)];
COG0 = COG(:,1);
COG2 = COG(:,2);
%%%uuuugh so NPSAD doesn't have the right cognitive measures. So I need to
%%%funpack extract. create file of the 4 cognitive measures
%%%cogIDs = [4282, 20016, 20023, 399];
%%%writematrix(cogIDs,'cogIDs.csv') 
%%% cogweighted = 0.716*NumericMemory + 0.774*Reasoning + -0.491*log(ReactionTime) + -0.495*log(VisualMemory);
%%now gotta run anova on cog that's been grouped by clinical group
%%
[~, ia] = intersect(COG(:,1), anhedonia_eid(:,2)); cogs{1} = COG(ia, 2);
[~, im] = intersect(COG(:,1), mood_eid(:,2)); cogs{2} = COG(im, 2);
[~, is] = intersect(COG(:,1), somatic_eid(:,2)); cogs{3} = COG(is, 2);
[~, ic] = intersect(COG(:,1), chronic_eid(:,2)); cogs{4} = COG(ic, 2);
[~, il] = intersect(COG(:,1), lateonset_eid(:,2)); cogs{5} = COG(il, 2);
[~, ise] = intersect(COG(:,1), severe_eid(:,2)); cogs{6} = COG(ise, 2);
[~, ih, jh] = intersect(COG(:,1), hetero_eid(:,2)); cogs{7} = COG(ih, 2);
for n = 1:6
    nans{n} = isnan(cogs{n});
end
for n = 1:6
    cog_nanless{n} = cogs{n}(~any(nans{n},2),:);
    kms_cog{n} = kms_TNM_3_all{n}(~any(nans{n},2),:);
end
%run the anova
for i=1:6
        p_cogs_6(i) = anova1(cog_nanless{i}, kms_cog{i},  'display', 'off');
end
%acute severity is statistically significant
%% If you had all 4 cog measures (you don't)
NumericMemory = [D.x4282_0_0 D.x4282_2_0];
Reasoning = [D.x20016_0_0 D.x20016_2_0];
ReactionTime = [D.x20023_0_0 D.x20023_2_0];
VisualMemory = [D.x399_0_0 D.x399_2_0];
cogweighted = 0.716*NumericMemory + 0.774*Reasoning + -0.491*log(ReactionTime) + -0.495*log(VisualMemory);
%okay that worked so change to the right codes and write it out correctly
%then remove nans and then you're good to run anova
%crap I don't know if I have the right codes
%% c-reactive protein
%load c-reactive protein (ID = 30710)
D = readtable('Subgroups.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);

[~, ia] = intersect(D.eid, anhedonia_eid(:,2)); B{1} = D.x30710_0_0(ia, :);
[~, im] = intersect(D.eid, mood_eid(:,2)); B{2} = D.x30710_0_0(im, :);
[~, is] = intersect(D.eid, somatic_eid(:,2)); B{3} = D.x30710_0_0(is, :);
[~, ic] = intersect(D.eid, chronic_eid(:,2)); B{4} = D.x30710_0_0(ic, :);
[~, il] = intersect(D.eid, lateonset_eid(:,2)); B{5} = D.x30710_0_0(il, :);
[~, ise] = intersect(D.eid, severe_eid(:,2)); B{6} = D.x30710_0_0(ise, :);
%remove missing data
for n = 1:6
    nans{n} = isnan(B{n});
end
for n = 1:6
    B_nanless{n} = B{n}(~any(nans{n},2),:);
    kms_B{n} = kms_TNM_3_all{n}(~any(nans{n},2),:);
end
%run anova
for i=1:6
        p_B_6(i) = anova1(B_nanless{i}, kms_B{i},  'display', 'off');
end
%none of them were significantly different
%% vascular
%get framingham 

%% characterize (no sig testing) on townsend index?


%% histograms for sample size
%%%%C-reactive Protein
figure
t = tiledlayout(3,2)
title(t, 'C-reactive protein')

nexttile
histogram(kms_B{1})
title('Anhedonia')

nexttile
histogram(kms_B{2})
title('Depressed Mood')

nexttile
histogram(kms_B{3})
title('Somatic Disturbance')

nexttile
histogram(kms_B{4})
title('Chronic')

nexttile
histogram(kms_B{5})
title('Late Onset')

nexttile
histogram(kms_B{6})
title('Acute Severity')

%%%% Cognition
figure
t = tiledlayout(3,2)
title(t, 'Cognition')

nexttile
histogram(kms_cog{1})
title('Anhedonia')

nexttile
histogram(kms_cog{2})
title('Depressed Mood')

nexttile
histogram(kms_cog{3})
title('Somatic Disturbance')

nexttile
histogram(kms_cog{4})
title('Chronic')

nexttile
histogram(kms_cog{5})
title('Late Onset')

nexttile
histogram(kms_cog{6})
title('Acute Severity')


%%%%% PRS_AD
figure
t = tiledlayout(3,2)
title(t, 'PRS Alzheimer')

nexttile
histogram(kms_PRSAD{1})
title('Anhedonia')

nexttile
histogram(kms_PRSAD{2})
title('Depressed Mood')

nexttile
histogram(kms_PRSAD{3})
title('Somatic Disturbance')

nexttile
histogram(kms_PRSAD{4})
title('Chronic')

nexttile
histogram(kms_PRSAD{5})
title('Late Onset')

nexttile
histogram(kms_PRSAD{6})
title('Acute Severity')

%%% PRS MDD
figure
t = tiledlayout(3,2)
title(t, 'PRS MDD')

nexttile
histogram(kms_PRSMDD{1})
title('Anhedonia')

nexttile
histogram(kms_PRSMDD{2})
title('Depressed Mood')

nexttile
histogram(kms_PRSMDD{3})
title('Somatic Disturbance')

nexttile
histogram(kms_PRSMDD{4})
title('Chronic')

nexttile
histogram(kms_PRSMDD{5})
title('Late Onset')

nexttile
histogram(kms_PRSMDD{6})
title('Acute Severity')

%%% Framingham
figure
t = tiledlayout(3,2)
title(t, 'Framingham')

nexttile
histogram(kms_fram{1})
title('Anhedonia')

nexttile
histogram(kms_fram{2})
title('Depressed Mood')

nexttile
histogram(kms_fram{3})
title('Somatic Disturbance')

nexttile
histogram(kms_fram{4})
title('Chronic')

nexttile
histogram(kms_fram{5})
title('Late Onset')

nexttile
histogram(kms_fram{6})
title('Acute Severity')