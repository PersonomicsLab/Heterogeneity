load('5aug22workspace.mat')

for i = 1:6
    for x = 1:length(kms_final{i})
        if kms_final{i}(x) == 1
        color{i}(x,:) = [0 0.4470 0.7410];
        else
        color{i}(x,:) = [0.8500 0.3250 0.0980];
        end
    end
end 

%% visualize PCA for all groups
figure
t = tiledlayout(2,3);
title(t, 'PCA for Clusters')
nexttile
scatter(score_final{1}(:,1), score_final{1}(:,2), [], color{1})
% set(gca,'XLim',[-30 30],'XTick',[-30:10:30], 'YLim',[-30 30],'YTick',[-30:10:30])
title('A. Anhedonia')
annotation('textbox',[.1 .78 .1 .1],'String','ARI = ','FitBoxToText','on');
xlabel('PCA Component 1')
ylabel('PCA Component 2')

nexttile
scatter(score_final{2}(:,1), score_final{2}(:,2), [], color{2})
% set(gca,'XLim',[-30 30],'XTick',[-30:10:30], 'YLim',[-30 30],'YTick',[-30:10:30])
title('B. Depressed Mood')
annotation('textbox',[.4 .78 .1 .1],'String','ARI = ','FitBoxToText','on');
xlabel('PCA Component 1')
ylabel('PCA Component 2')

nexttile
scatter(score_final{3}(:,1), score_final{3}(:,2), [], color{3})
% set(gca,'XLim',[-30 30],'XTick',[-30:10:30], 'YLim',[-30 30],'YTick',[-30:10:30])
title('C. Somatic Disturbance')
annotation('textbox',[.7 .78 .1 .1],'String','ARI = ','FitBoxToText','on');
xlabel('PCA Component 1')
ylabel('PCA Component 2')

nexttile
scatter(score_final{4}(:,1), score_final{4}(:,2), [], color{4})
% set(gca,'XLim',[-30 30],'XTick',[-30:10:30], 'YLim',[-30 30],'YTick',[-30:10:30])
title('D. Chronic')
annotation('textbox',[.1 .3 .1 .1],'String','ARI = , sig. Cognition','FitBoxToText','on');
xlabel('PCA Component 1')
ylabel('PCA Component 2')

nexttile
scatter(score_final{5}(:,1), score_final{5}(:,2), [], color{5})
% set(gca,'XLim',[-30 30],'XTick',[-30:10:30], 'YLim',[-30 30],'YTick',[-30:10:30])
title('E. Late Onset')
annotation('textbox',[.4 .3 .1 .1],'String','ARI = ','FitBoxToText','on');
xlabel('PCA Component 1')
ylabel('PCA Component 2')

nexttile
scatter(score_final{6}(:,1), score_final{6}(:,2), [], color{6})
% set(gca,'XLim',[-30 30],'XTick',[-30:10:30], 'YLim',[-30 30],'YTick',[-30:10:30])
title('F. Acute Severe')
annotation('textbox',[.7 .3 .1 .1],'String','ARI = , sig. Cognition','FitBoxToText','on');
xlabel('PCA Component 1')
ylabel('PCA Component 2')


%% ARI error bars
for i = 1:6
    temp = cell2mat(ARI_all(i, :));
    temp(temp==0)=NaN;
    ARISEM{i} = nanstd(temp, [], 'all')/sqrt(900);
end

ARImean = cell2mat(ARI_dist);
ARIsem = cell2mat(ARISEM);
errlow = ARImean - ARIsem;
errhigh = ARImean + ARIsem;
%%
scoretotal = zeros(length(idxTrain{2,1, 6}),31);
scoretotal(idxTrain{2,1,6}, :) = score_temp{2,1, 6}(:,1:31);
scoretotal(idxNew{2,1,6}, :) = score_new{2,1, 6}; 
kmstotal = zeros(length(idxTrain{2,1,6}),1);
kmstotal(idxTrain{2,1,6}, :) = kms_temp{2,1,6};
kmstotal(idxNew{2,1,6}, :) = cl_idx{2,1,6};
scoretotal2 = zeros(length(idxTrain{4,1, 6}),31);
scoretotal2(idxTrain{4,1,6}, :) = score_temp{4,1, 6}(:,1:31);
scoretotal2(idxNew{4,1,6}, :) = score_new{4,1, 6}; 
    for x = 1:length(kmstotal)
        if kmstotal(x) == 1
        color(x,:) = [0.4824    0.6902    0.7804];
        else
        color(x,:) = [0.5686    0.5569    0.5569];
        end
    end
%% aim 2 figure draft 1
% load('PCA_ARI_allage_75jul8.mat')
figure
t = tiledlayout(2,2);

nexttile % A. example PCA
% plot one of the folds 
scatter(scoretotal(:,1), scoretotal(:,2),[], color)
xlabel('PCA Component 1')
ylabel('PCA Component 2')
t = annotation('textbox', [0.38, 0.78, 0.1, 0.1], 'String', ['Cluster 1' newline]);
t.Color = [0.5686    0.5569    0.5569];
t = annotation('textbox', [0.38, 0.78, 0.1, 0.1], 'String', [newline 'Cluster 2']);
t.Color = [0.4824    0.6902    0.7804];
title('A. Example Clustering for Acute Impairment')
set(gca,'FontSize',14)

nexttile % B. cluster stability 
% plot one of the other folds
scatter(scoretotal2(:,1), scoretotal2(:,2),[], color)
t = annotation('textbox', [0.84, 0.78, 0.1, 0.1], 'String', ['Cluster 1' newline]);
t.Color = [0.5686    0.5569    0.5569];
t = annotation('textbox', [0.84, 0.78, 0.1, 0.1], 'String', [newline 'Cluster 2']);
t.Color = [0.4824    0.6902    0.7804];
xlabel('PCA Component 1')
ylabel('PCA Component 2')
title('B. Example of Cluster Stability for Acute Impairment')
set(gca,'FontSize',14)
%%%okay I figured out the issue - the order for the subjects is different
%%%for every iteration

nexttile % C. ARI of all groups
% bar plot ARI_dist {0.712336147447252,0.405592982808764,0.753152660997236,0.611499227608624,0.885657983708292,0.868676573529961}
% ARInew = [ARImean(1:4), ARImean(6), ARImean(5)];
% ARIsemnew = [ARIsem(1:4), ARIsem(6), ARIsem(5)];
bar([1:6], ARInew)
text([1:length(ARInew)], ARInew', num2str(ARInew','%0.2f'),'HorizontalAlignment','center','VerticalAlignment','bottom')
set(gca,'FontSize',14)
ylabel('Adjusted Rand Index')
ylim([0, 1])
xticklabels({'Anhedonia', 'Depressed Mood', 'Somatic Disturbance', 'Chronicity', 'Acute Impairment', 'Late Onset'})
hold on
er = errorbar([1:6],ARInew, ARIsemnew,ARIsemnew);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
title('C. Cluster Stability For All Clinically Homogeneous Groups')
hline(.03)
t = annotation('textbox', [0.07, 0.088, 0.1, 0.1], 'String', "Null Threshold", 'EdgeColor', 'none', 'FontSize',12);
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
title('D. Cognition Cluster Validation For All Clinically Homogeneous Groups')
set(gca,'FontSize',14)
t = annotation('textbox', [0.78, 0.3, 0.1, 0.1], 'String', "**");
t.LineStyle = 'none';
t.FontSize = 24;
t = annotation('textbox', [0.84, 0.3, 0.1, 0.1], 'String', "**");
t.LineStyle = 'none';
t.FontSize = 24;