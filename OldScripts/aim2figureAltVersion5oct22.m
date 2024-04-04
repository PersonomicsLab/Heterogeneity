figure
t = tiledlayout(1,2);

nexttile % C. ARI of all groups
% bar plot ARI_dist {0.712336147447252,0.405592982808764,0.753152660997236,0.611499227608624,0.885657983708292,0.868676573529961}
% ARInew = [ARImean(1:4), ARImean(6), ARImean(5)];
% ARIsemnew = [ARIsem(1:4), ARIsem(6), ARIsem(5)];
bar([1:6], ARInew)
text([1:length(ARInew)], ARInew', num2str(ARInew','%0.2f'),'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontSize', 14)
set(gca,'FontSize',14)
ylabel('Adjusted Rand Index')
ylim([0, 1])
xticklabels({'Anhedonia', 'Depressed Mood', 'Somatic Disturbance', 'Chronicity', 'Acute Impairment', 'Late Onset'})
hold on
er = errorbar([1:6],ARInew, ARIsemnew,ARIsemnew);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
title('Cluster Stability')
hline(.03)
t = annotation('textbox', [0.07, 0.15, 0.1, 0.1], 'String', "Null Threshold", 'EdgeColor', 'none', 'FontSize',12);
nexttile % D. Cognition violin plots fo clusters for all groups
% find the cognition data wherever that is and find the violin graph code
% wherever that is
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions/violin (1)')
space = 3;
pos1 = [2 5 8 11 14 17];
pos2 = [2.7 5.7 8.7 11.7 14.7 17.7];

hold on
violin({cog_nanless{1}(kms_cog{1}==1, :), cog_nanless{2}(kms_cog{2}==1, :), cog_nanless{3}(kms_cog{3}==1, :), cog_nanless{4}(kms_cog{4}==1, :), cog_nanless{6}(kms_cog{6}==1, :), cog_nanless{5}(kms_cog{5}==1, :)}, 'x', pos1 ,'facecolor', [0.5 0.5 0.5], 'medc', '')
violin({cog_nanless{1}(kms_cog{1}==2, :), cog_nanless{2}(kms_cog{2}==2, :), cog_nanless{3}(kms_cog{3}==2, :), cog_nanless{4}(kms_cog{4}==2, :), cog_nanless{6}(kms_cog{6}==2, :), cog_nanless{5}(kms_cog{5}==2, :)}, 'x', pos2 ,'facecolor', [0.6 0.8 1], 'medc', '')

hold off
set(gca,'xticklabel',{'Anhedonia', 'Depressed Mood', 'Somatic Disturbance', 'Chronicity', 'Acute Impairment', 'Late Onset'})
xpos = mean([pos1; pos2],1);
xticks(xpos)
xlim([0, max(xpos)+space/2])
ylabel('Cognition Score')
title('Cognition Cluster Validation')
set(gca,'FontSize',14)
t = annotation('textbox', [0.78, 0.75, 0.1, 0.1], 'String', "**");
t.LineStyle = 'none';
t.FontSize = 24;
t = annotation('textbox', [0.85, 0.75, 0.1, 0.1], 'String', "**");
t.LineStyle = 'none';
t.FontSize = 24;