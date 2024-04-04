%% find age, sex, and townsend then copy the graph code
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data/allage')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Code')

%%load clinical group eids
anhedonia_eid = table2array(readtable('HomogeneousGroups_anhedonia.csv'));
mood_eid = table2array(readtable('HomogeneousGroups_mood.csv'));
somatic_eid = table2array(readtable('HomogeneousGroups_somatic.csv'));
chronic_eid = table2array(readtable('HomogeneousGroups_chronic.csv'));
lateonset_eid = table2array(readtable('HomogeneousGroups_lateonset.csv'));
severe_eid = table2array(readtable('HomogeneousGroups_severe.csv'));
hetero_eid = table2array(readtable('HomogeneousGroups_HeterogeneousComparison.csv'));
% anhdep_eid = table2array(readtable('HomogeneousGroups_anhdep.csv'));

%% age
%where is age
D = readtable('Subgroups_new.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);

[~, ia] = intersect(D.eid, anhedonia_eid(:,1)); A{1} = D.x21003_2_0(ia, :);
[~, im] = intersect(D.eid, mood_eid(:,1)); A{2} = D.x21003_2_0(im, :);
[~, is] = intersect(D.eid, somatic_eid(:,1)); A{3} = D.x21003_2_0(is, :);
[~, ic] = intersect(D.eid, chronic_eid(:,1)); A{4} = D.x21003_2_0(ic, :);
[~, il] = intersect(D.eid, lateonset_eid(:,1)); A{5} = D.x21003_2_0(il, :);
[~, ise] = intersect(D.eid, severe_eid(:,1)); A{6} = D.x21003_2_0(ise, :);
%remove missing data
for n = 1:6
    nans{n} = isnan(A{n});
end
for n = 1:6
    A_nanless{n} = A{n}(~any(nans{n},2),:);
    kms_A{n} = kms_final{n}(~any(nans{n},2),:);
end
%run anova
for i=1:6
    p_A(i) = anova1(A_nanless{i}, kms_A{i},  'display', 'off');
end
[pthr,pcor,padj] = fdr(p_A);
sigA = find(pcor < 0.0500);

%% sex
%where is sex
D = readtable('kayla.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);

[~, ia] = intersect(D.eid, anhedonia_eid(:,1)); S{1} = D.x31_0_0(ia, :);
[~, im] = intersect(D.eid, mood_eid(:,1)); S{2} = D.x31_0_0(im, :);
[~, is] = intersect(D.eid, somatic_eid(:,1)); S{3} = D.x31_0_0(is, :);
[~, ic] = intersect(D.eid, chronic_eid(:,1)); S{4} = D.x31_0_0(ic, :);
[~, il] = intersect(D.eid, lateonset_eid(:,1)); S{5} = D.x31_0_0(il, :);
[~, ise] = intersect(D.eid, severe_eid(:,1)); S{6} = D.x31_0_0(ise, :);
%remove missing data
for n = 1:6
    nans{n} = isnan(S{n});
end
for n = 1:6
    S_nanless{n} = S{n}(~any(nans{n},2),:);
    kms_S{n} = kms_final{n}(~any(nans{n},2),:);
end
for n = 1:6
[h_S(n),p_S(n)]=chi2gof(S_nanless{n}(kms_S{n}));  %%result=1,3,5
end
% for n = 1:6
% [h_S(n),p_S(n)]=ttest2((S_nanless{n}(kms_S{n})==1), (S_nanless{n}(kms_S{n})==1)); %%result=1,3,4,5,6
% end

% for n = 1:6
%     [h_S(n),p_S(n)]=fishertest([scounta(n,:); scountb(n,:)]); %%result=1,3,4,5,6
% end
[pthr,pcor,padj] = fdr(p_S);
sigS = find(pcor < 0.0500);

%% townsend (dang, no cigar)
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
    kms_T{n} = kms_final{n}(~any(nans{n},2),:);
end
%run anova
for i=1:6
    p_T(i) = anova1(T_nanless{i}, kms_T{i},  'display', 'off');
end
[pthr,pcor,padj] = fdr(p_T);
sigT = find(pcor < 0.0500);


%% cognition
D = readtable('Cognition.tsv','FileType','text');
D = standardizeMissing(D,-3); D = standardizeMissing(D,-1); D = standardizeMissing(D,-818); D = standardizeMissing(D,-121); D = standardizeMissing(D,-7);
Reasoning = [D.x20016_0_0 D.x20016_2_0];
ReactionTime = [D.x20023_0_0 D.x20023_2_0];
COG = [D.eid 0.774*Reasoning + -0.491*log(ReactionTime)];
COG0 = COG(:,1);
COG2 = COG(:,2);

[~, ia] = intersect(COG(:,1), anhedonia_eid(:,1)); cogs{1} = COG(ia, 2);
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
    p_cogs_6(i) = anova1(cog_nanless{i}, kms_cog{i},  'display', 'off');
end
[pthr,pcor,padj] = fdr(p_cogs_6);
sigC = find(pcor < 0.0500);

%% visualize template
figure
t = tiledlayout(2,2);

nexttile
tit = 'Cognition';
label = 'Cognition Score';
type = 'violin';
visualize(type, cog_nanless, kms_cog, tit, label)

nexttile
tit = 'Townsend';
label = 'Townsend Score';
type = 'violin';
visualize(type, T_nanless, kms_T, tit, label)

nexttile
tit = 'Age';
label = 'Age (years)';
type = 'violin';
visualize(type, A_nanless, kms_A, tit, label)

nexttile
tit = 'Sex';
label = 'Sex (%)';
type = 'bar';
visualize(type, S_nanless, kms_S, tit, label)

%% Functions
for i = 1:6
    [zstat(i)] = myztest(scounta(i,1),scountb(i,1),(scounta(i,1)+scounta(i,2)),(scountb(i,1)+scountb(i,2)));
end
function visualize(type, phen_nanless, phen_kms, tit, label)
if type == "bar"
    for i = 1:6
        counta(i,1) = sum(phen_nanless{i}(phen_kms{i}==1, :)==0)/length(phen_nanless{i}(phen_kms{i}==1, :))*100;
        counta(i,2) = sum(phen_nanless{i}(phen_kms{i}==1, :)==1)/length(phen_nanless{i}(phen_kms{i}==1, :))*100;
        countb(i,1) = sum(phen_nanless{i}(phen_kms{i}==2, :)==0)/length(phen_nanless{i}(phen_kms{i}==2, :))*100;
        countb(i,2) = sum(phen_nanless{i}(phen_kms{i}==2, :)==1)/length(phen_nanless{i}(phen_kms{i}==2, :))*100;
    end
    space = 5;
    pos1 = [1.5 6 10.5 15 19.5 24];
    pos2 = [3 7.5 12 16.5 21 25.5];
    % space = 3;
    % pos1 = [2 5 8 11 14 17];
    % pos2 = [2.7 5.7 8.7 11.7 14.7 17.7];
    
    hold on
    ba = bar(pos1, counta, 0.3, 'stacked')
    set(ba,{'FaceColor'},{[0.6510    0.6510    0.6510]; [0.9412    0.9412    0.9412]});
    ba = bar(pos2, countb, 0.3, 'stacked')
    set(ba,{'FaceColor'},{[0.6941    0.8471    0.9804];[0.8431    0.9176    0.9804]});
    hold off
    
    legend('Cluster 1 Female', 'Cluster 1 Male', 'Cluster 2 Female', 'Cluster 2 Male', 'Location', 'northwest')
    set(gca,'xticklabel',{'Anhedonia', 'Depressed Mood', 'Somatic Disturbance', 'Chronicity', 'Acute Impairment', 'Late Onset'})
    xpos = mean([pos1; pos2],1);
    xticks(xpos)
    xlim([0, max(xpos)+space/2])
    ylabel(label)
    title(tit)
    set(gca,'FontSize',14)
    hline(50)
    
else
    addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions/violin (1)')
    space = 3;
    pos1 = [2 5 8 11 14 17];
    pos2 = [2.7 5.7 8.7 11.7 14.7 17.7];
    
    hold on
    violin({phen_nanless{1}(phen_kms{1}==1, :), phen_nanless{2}(phen_kms{2}==1, :), phen_nanless{3}(phen_kms{3}==1, :), phen_nanless{4}(phen_kms{4}==1, :), phen_nanless{6}(phen_kms{6}==1, :), phen_nanless{5}(phen_kms{5}==1, :)}, 'x', pos1 ,'facecolor', [0.5 0.5 0.5], 'medc', '')
    violin({phen_nanless{1}(phen_kms{1}==2, :), phen_nanless{2}(phen_kms{2}==2, :), phen_nanless{3}(phen_kms{3}==2, :), phen_nanless{4}(phen_kms{4}==2, :), phen_nanless{6}(phen_kms{6}==2, :), phen_nanless{5}(phen_kms{5}==2, :)}, 'x', pos2 ,'facecolor', [0.6 0.8 1], 'medc', '')
    hold off
    
    set(gca,'xticklabel',{'Anhedonia', 'Depressed Mood', 'Somatic Disturbance', 'Chronicity', 'Acute Impairment', 'Late Onset'})
    xpos = mean([pos1; pos2],1);
    xticks(xpos)
    xlim([0, max(xpos)+space/2])
    ylabel(label)
    title(tit)
    set(gca,'FontSize',14)
end
% t = annotation('textbox', [0.78, 0.3, 0.1, 0.1], 'String', "**");
% t.LineStyle = 'none';
% t.FontSize = 24;
% t = annotation('textbox', [0.84, 0.3, 0.1, 0.1], 'String', "**");
% t.LineStyle = 'none';
% t.FontSize = 24;
end

function [zstat] = myztest(x1,x2,n1,n2)
%Z test for the difference in two proportions
%   Z test to compare two different proportions, using the x1 and x2 as the
%   successes of each group and n1 and n2 as the total in the group. Answer
%   = the Z stat, in which the null hypthothesis is rejected if Zstat<-1.96
%   or if Z>1.96 when alpha = 0.05
x=((x1/n1)-(x2/n2))-(0);
y=sqrt((((x1+x2)/(n1+n2))*(1-((x1+x2)/(n1+n2))))*((1/n1+1/n2)));
zstat=x/y;
end