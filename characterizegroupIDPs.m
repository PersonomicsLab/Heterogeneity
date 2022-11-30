%% run first two sections of PCAforallage.m

%%
for i = 1:7
 meanZ(:,i) = mean(Z{i},1);
end
%%figure for all FA of all groups (1-19)
%  %%%% calculate mean of each IDP for each clinical group (88 by 7 matrix)
% figure; 
% h = bar(meanZ(1:19, :));
% legend({'anhedonia','mood','somatic', 'chronic','late onset', 'severe', 'heterogenous'},'Location','EastOutside')
% ylabel('normative mean')
% xticks(1:19)
% xticklabels(IDPnames(1:19, 1))
% title('FAs')
% 
% %% test
% figure; 
% h = bar(aaa(:, 1:19));
% xticklabels({'anhedonia','mood','somatic', 'chronic','late onset', 'severe', 'heterogenous'},'Location','EastOutside')
% ylabel('normative mean')
% xticks(1:6)
% legend(IDPnames(1:19, 1))
% title('FAs')
% 
% %% figure for WMH (20)
% figure; 
% h = bar(meanZ(20, :));
% legend({'anhedonia','mood','somatic', 'chronic','late onset', 'severe', 'heterogenous'},'Location','EastOutside')
% ylabel('normative mean')
% xticks(1)
% xticklabels(IDPnames(20, 1))
% title('WMH')
% %% GMV (21-48)
% figure; 
% h = bar(meanZ(21:48, :));
% legend({'anhedonia','mood','somatic', 'chronic','late onset', 'severe', 'heterogenous'},'Location','EastOutside')
% ylabel('normative mean')
% xticks(1:28)
% xticklabels(IDPnames(21:48, 1))
% title('GMVs')
% 
% %% CT (49-60)
% figure; 
% h = bar(meanZ(49:60, :));
% legend({'anhedonia','mood','somatic', 'chronic','late onset', 'severe', 'heterogenous'},'Location','EastOutside')
% ylabel('normative mean')
% xticks(1:12)
% xticklabels(IDPnames(49:60, 1))
% title('CTs')
% 
% %% FC (61-88)
% figure; 
% h = bar(meanZ(61:88, :));
% legend({'anhedonia','mood','somatic', 'chronic','late onset', 'severe', 'heterogenous'},'Location','EastOutside')
% ylabel('normative mean')
% xticks(1:28)
% xticklabels(IDPnames(61:88, 1))
% title('FCs')

%%Consolidated IDPs
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions/github_repo (1)')
%seperate by modality
FA = [mean(meanZ(1:3, :),1);  mean(meanZ(4:5, :),1); mean(meanZ(6:11, :),1); mean(meanZ(12:13, :),1); mean(meanZ(14:15, :),1); mean(meanZ(16:17, :),1); mean(meanZ(18:19, :),1)];
FAname = [{'Corpus Collosum', 'ALIC', 'Corona Radiata', 'Cingulum Cingulate Gyrus', 'Cingulum Hippocampus', 'Superior Longitudinal Fasciculus', 'Uncinate Fasciculus'}];

GMV = [mean(meanZ(21:22, :),1); mean(meanZ(23:26, :),1); mean(meanZ(27:32, :),1);  mean(meanZ(33:34, :),1); mean(meanZ(35:38, :),1);  mean(meanZ(39:40, :),1);  mean(meanZ(41:42, :),1); mean(meanZ(43:44, :),1); mean(meanZ(45:46, :),1); mean(meanZ(47:48, :),1)];
GMVname = [{'Insular Cortex', 'Superior Temporal Gyrus', 'Middle Temporal Gyrus', 'Frontal Medial Cortex', 'Cingulate Gyrus', 'Frontal Orbital Cortex', 'Thalamus', 'Hippocampus', 'Amygdala', 'Ventral Striatum'}];

CT = [mean(meanZ([49,55], :),1);  mean(meanZ([50,56], :),1); mean(meanZ([51,57], :),1); mean(meanZ([52,58], :),1); mean(meanZ([53, 59], :),1)];
CTname = [{'Entorhinal', 'Fusiform', 'Inferiorparietal', 'Parahippocampal', 'Superior Frontal'}];

FC = [mean(meanZ(61:63, :),1); mean(meanZ(64:65, :),1); mean(meanZ(66:67, :),1); mean(meanZ([68,69,74], :),1); mean(meanZ([70,71,75,76,79,80], :),1); mean(meanZ([72,73,77,78,81,82], :),1); meanZ(83, :); mean(meanZ([84,85,86,87], :),1); meanZ(88, :)];
FCname = [{'DMN Amplitude', 'CEN Amplitude', 'Salience Network Amplitude', 'DMN/DMN Connectivity', 'DMN/CEN Connectivity', 'DMN/SN Connectivity', 'CEN/CEN Connectivity', 'CEN/SN Connectivity', 'SN/SN Connectivity'}];

%%sig IDPS
sigIDPs = [{'Anterior Corona Radiata FA', 'Superior Longitudinal Fasciculus FA', 'White Matter Hyperintensity', 'Insular Cortex GMV', 'Superior Temporal Gyrus GMV', 'Cingulate Gyrus GMV', 'Frontal Orbital Cortex GMV', 'Amygdala GMV', 'Entorhinal CT', 'Fusiform CT', 'Inferiorparietal CT', 'Superiorparietal CT', 'DMN/DMN Connectivity', 'CEN/SN Connectivity'}];
sigmeanZ = [meanZ(sigidx([1:5]), :); mean(meanZ(sigidx([6,7]), :),1); meanZ(sigidx(8), :); mean(meanZ(sigidx([9,10]), :),1); meanZ(sigidx([11:end]), :)];
%%SEM = STD(group)/sqrt(lengthofgroup)
for i = 1:7
    for x = 1:88
        SEM(x,i) = std(Z{i}(:,x))/sqrt(length(Z{i}));
    end
end
SEMkeep = [SEM(sigidx([1:5]), :); mean(SEM(sigidx([6,7]), :),1); SEM(sigidx(8), :); mean(SEM(sigidx([9,10]), :),1); SEM(sigidx([11:end]), :)];

%% Make spider plots
% figure
% t = tiledlayout(2,4);
% title(t, 'Imaging Features for each Clinical Group')
% % legend(t, {'anhedonia','mood','somatic', 'chronic','late onset', 'severe'})
% nexttile 
% spider_plot(FA(:,1:6)','AxesLabels', FAname, 'LabelFontSize', 8, 'AxesLimits', [-.2, -.2, -.2, -.2, -.2, -.2, -.2, -.1; .1, .1, .1, .1, .1, .1, .1, .3]);
% title('A. White Matter') 
% nexttile
% spider_plot(GMV(:,1:6)', 'AxesLabels', GMVname, 'LabelFontSize', 8, 'AxesLimits', [-.15, -.15, -.15, -.15, -.15, -.15, -.15, -.15, -.15, -.15; .15, .15, .15, .15, .15, .15, .15, .15, .15, .15]);
% title('B. Gray Matter Volume')
% nexttile
% spider_plot(CT(:,1:6)', 'AxesLabels', CTname, 'LabelFontSize', 8, 'AxesLimits', [-.15, -.15, -.15, -.15, -.15; .1, .1, .1, .1, .1])
% title('C. Cortical Thickness')
% nexttile
% spider_plot(FC(:,1:6)', 'AxesLabels', FCname, 'LabelFontSize', 8, 'AxesLimits', [-.1, -.1, -.1, -.1, -.1, -.1, -.1, -.1, -.1; .15, .15, .15, .15, .15, .15, .15, .15, .15])
% title('D. Resting State fMRI')
% 
% nexttile([1 4])
% bar(sigmeanZ);
% ylabel('Normative Deviation Mean')
% set(gca,'XLim',[0 15],'XTick',[1:1:14])
% xticklabels(sigIDPs)
% vline([1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5]);
% title('E. Imaging Measures Significantly Different between Clinical Groups')
% hold on
% % Find the number of groups and the number of bars in each group
% [ngroups, nbars] = size(sigmeanZ);
% % Calculate the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% % Set the position of each error bar in the centre of the main bar
% % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
% for i = 1:nbars
%     % Calculate center of each bar
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     er = errorbar(x, sigmeanZ(:,i), SEMkeep(:,i), 'k', 'linestyle', 'none');
% end
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 
% hold off
% legend({'Anhedonia','Low Mood','Somatic', 'Chronic','Late Onset', 'Acute Severe', 'Heterogenous'},'Location','EastOutside')

%% different (consistent) axes
figure
% t = tiledlayout(2,4);
% % title(t, 'Imaging Features for each Clinical Group')
% nexttile 
% spider_plot(FC(:,1:6)', 'AxesLabels', FCname, 'LabelFontSize',13, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
% title('A. Resting State fMRI')
% set(gca,'FontSize',15)
% nexttile
% spider_plot(GMV(:,1:6)', 'AxesLabels', GMVname, 'LabelFontSize',13, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none');
% title('B. Gray Matter Volume')
% set(gca,'FontSize',15)
% nexttile
% spider_plot(CT(:,1:6)', 'AxesLabels', CTname, 'LabelFontSize', 13, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
% title('C. Cortical Thickness')
% set(gca,'FontSize',15)
% nexttile
% spider_plot(FA(:,1:6)','AxesLabels', FAname, 'LabelFontSize', 13, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none');
% title('D. Fractional Anisotropy') 
% set(gca,'FontSize',15)
% % lg  = legend({'anhedonia','mood','somatic', 'chronic','late onset', 'severe'}); 
% % lg.Layout.Tile = 'South';
% 
% nexttile([1 4])
bar(sigmeanZ);
ylabel('Normative Deviation Mean')
ylim([-.3 .31])
set(gca,'XLim',[0 15],'XTick',[1:1:14])
xticklabels(sigIDPs)
set(gca,'FontSize',18)
%t = title('E. Imaging Measures Significantly Different between Clinically Isolated Groups')
%t.FontSize = 22;
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(sigmeanZ);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, sigmeanZ(:,i), SEMkeep(:,i), 'k', 'linestyle', 'none');
end
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
vline([1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5]);
legend({'Anhedonia','Low Mood','Somatic', 'Chronic','Late Onset', 'Acute Impairment', 'Heterogenous'},'Location','Northeast')