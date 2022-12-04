%%%run first two sections of PCAforallage.m

%% Prep data and labels
for i = 1:7
 meanZ(:,i) = mean(Z{i},1);
end

%%%Consolidated IDPs
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions/github_repo (1)')
%seperate by modality
FA = [mean(meanZ(1:3, :),1);  mean(meanZ(4:5, :),1); mean(meanZ(6:11, :),1); mean(meanZ(12:13, :),1); mean(meanZ(14:15, :),1); mean(meanZ(16:17, :),1); mean(meanZ(18:19, :),1)];
FAname = [{'Corpus Collosum', 'ALIC', 'Corona Radiata', 'Cingulum Cingulate Gyrus', 'Cingulum Hippocampus', 'Superior Longitudinal Fasciculus', 'Uncinate Fasciculus'}];

GMV = [mean(meanZ(21:22, :),1); mean(meanZ(23:26, :),1); mean(meanZ(27:32, :),1);  mean(meanZ(33:34, :),1); mean(meanZ(35:38, :),1);  mean(meanZ(39:40, :),1);  mean(meanZ(41:42, :),1); mean(meanZ(43:44, :),1); mean(meanZ(45:46, :),1); mean(meanZ(47:48, :),1)];
GMVname = [{'Insular Cortex', 'Superior Temporal Gyrus', 'Middle Temporal Gyrus', 'Frontal Medial Cortex', 'Cingulate Gyrus', 'Frontal Orbital Cortex', 'Thalamus', 'Hippocampus', 'Amygdala', 'Ventral Striatum'}];

CT = [mean(meanZ([49,55], :),1);  mean(meanZ([50,56], :),1); mean(meanZ([51,57], :),1); mean(meanZ([52,58], :),1); mean(meanZ([53, 59], :),1)];
CTname = [{'Entorhinal', 'Fusiform', 'Inferiorparietal', 'Parahippocampal', 'Superior Frontal'}];

FC = [mean(meanZ(61:63, :),1); mean(meanZ(64:65, :),1); mean(meanZ(66:67, :),1); mean(meanZ([68,69,74], :),1); mean(meanZ([70,71,75,76,79,80], :),1); mean(meanZ([72,73,77,78,81,82], :),1); meanZ(83, :); mean(meanZ([84,85,86,87], :),1); meanZ(88, :)];
FCname = [{'DMN Amplitude', 'CEN Amplitude', 'SN Amplitude', 'DMN/DMN Connectivity', 'DMN/CEN Connectivity', 'DMN/SN Connectivity', 'CEN/CEN Connectivity', 'CEN/SN Connectivity', 'SN/SN Connectivity'}];

%%sig IDPS
% sigIDPs = [{'Right Anterior Corona Radiata FA', 'Right Superior Longitudinal Fasciculus FA', 'White Matter Hyperintensity', 'Cingulate Gyrus GMV', 'Amygdala GMV', 'Insular Cortex GMV', 'Superior Temporal Gyrus GMV', 'Frontal Orbital Cortex GMV', 'Entorhinal CT', 'Fusiform CT', 'Inferiorparietal CT', 'Superiorparietal CT', 'DMN/DMN Connectivity', 'CEN/SN Connectivity'}];
sigIDPs = [{'Right Anterior Corona Radiata FA', 'Right Superior Longitudinal Fasciculus FA', 'White Matter Hyperintensity', 'Left Cingulate Gyrus GMV', 'Right Cingulate Gyrus GMV', 'Left Amygdala GMV', 'Right Amygdala GMV', 'Left Insular Cortex GMV', 'Left Superior Temporal Gyrus GMV', 'Left Frontal Orbital Cortex GMV', 'Left Entorhinal CT', 'Left Fusiform CT', 'Left Inferiorparietal CT', 'Left Superiorparietal CT', 'DMN/DMN Connectivity', 'CEN/SN Connectivity'}];
sigmeanZ = [meanZ(sigidx([1:3]), :); meanZ(sigidx([6,7]), :); meanZ(sigidx([9,10]), :); meanZ(sigidx([4,5]), :); meanZ(sigidx(8), :); meanZ(sigidx([11:end]), :)];
%%SEM = STD(group)/sqrt(lengthofgroup)
for i = 1:7
    for x = 1:88
        SEM(x,i) = std(Z{i}(:,x))/sqrt(length(Z{i}));
    end
end
sigSEM = [SEM(sigidx([1:3]), :); SEM(sigidx([6,7]), :); SEM(sigidx([9,10]), :); SEM(sigidx([4,5]), :); SEM(sigidx(8), :); SEM(sigidx([11:end]), :)];


%% fisualization, axes between modalities are consistent
figure;
t = tiledlayout(2,4);
% title(t, 'Imaging Features for each Clinical Group')
nexttile 
spider_plot(FC(:,1:6)', 'AxesLabels', FCname, 'LabelFontSize',12, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
%%% https://www.mathworks.com/matlabcentral/fileexchange/59561-spider_plot
title('A. Resting State fMRI')
set(gca,'FontSize',14)
nexttile
spider_plot(GMV(:,1:6)', 'AxesLabels', GMVname, 'LabelFontSize',12, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none');
title('B. Gray Matter Volume')
set(gca,'FontSize',14)
nexttile
spider_plot(CT(:,1:6)', 'AxesLabels', CTname, 'LabelFontSize', 12, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
title('C. Cortical Thickness')
set(gca,'FontSize',14)
nexttile
spider_plot(FA(:,1:6)','AxesLabels', FAname, 'LabelFontSize', 12, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none');
title('D. Fractional Anisotropy')
set(gca,'FontSize',14)
% lg  = legend({'anhedonia','mood','somatic', 'chronic','late onset', 'severe'}); 
% lg.Layout.Tile = 'South';

nexttile([1 4])
bar(sigmeanZ);
ylabel('Normative Deviation Mean')
ylim([-.3 .31])
% ax = axes('Parent',H);
set(gca,'XLim',[0 17],'XTick',1:1:16)
xticklabels(sigIDPs)
set(gca,'FontSize',16)
stats = stats_from_c(c);
disp(['found total of ', num2str(length(find(stats))),' significant pairs'])
[stats_Y,stats_X1,stats_X2] = plot_stats(gca,stats); %https://www.mathworks.com/matlabcentral/fileexchange/92613-statistical-significance
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
    er = errorbar(x, sigmeanZ(:,i), sigSEM(:,i), 'k', 'linestyle', 'none');
end
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
xline([2.5, 3.5, 10.5, 14.5], '--', 'LineWidth', 2)
hold off
vline([1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 15.5]);
%%%https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline?tab=reviews%2F2081773
legend({'Anhedonia','Low Mood','Somatic', 'Chronic','Late Onset', 'Acute Impairment', 'Heterogenous'},'Location','Northeast')


%% ended up saving them out individually for space reasons
figure
spider_plot(FC(:,1:6)', 'AxesLabels', FCname, 'LabelFontSize',14, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
%%% https://www.mathworks.com/matlabcentral/fileexchange/59561-spider_plot
title('A. Resting State fMRI')
set(gca,'FontSize',18)
set(gcf, 'Position', [550   557   715   420])
% print(gcf,'aim1figA_2dec22','-dpng','-r300');
exportgraphics(gcf,'transparentA.eps',...   % since R2020a
    'ContentType','vector',...
    'BackgroundColor','none')
close all

figure
spider_plot(GMV(:,1:6)', 'AxesLabels', GMVname, 'LabelFontSize',14, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none');
title('B. Gray Matter Volume')
set(gca,'FontSize',18)
set(gcf, 'Position', [550   557   715   420])
% print(gcf,'aim1figB_2dec22','-dpng','-r300');
exportgraphics(gcf,'transparentB.eps',...   % since R2020a
    'ContentType','vector',...
    'BackgroundColor','none')
close all

figure
spider_plot(CT(:,1:6)', 'AxesLabels', CTname, 'LabelFontSize', 14, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
title('C. Cortical Thickness')
set(gca,'FontSize',18)
set(gcf, 'Position', [550   557   715   420])
% print(gcf,'aim1figC_2dec22','-dpng','-r300');
exportgraphics(gcf,'transparentC.eps',...   % since R2020a
    'ContentType','vector',...
    'BackgroundColor','none')
close all

figure
spider_plot(FA(:,1:6)','AxesLabels', FAname, 'LabelFontSize', 14, 'AxesPrecision', 2, 'AxesLimits', [-.23, -.23, -.23, -.23, -.23, -.23, -.23; .12, .12, .12, .12, .12, .12, .12], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none');
title('D. Fractional Anisotropy')
set(gca,'FontSize',18)
set(gcf, 'Position', [550   557   715   420])
% print(gcf,'aim1figD_2dec22','-dpng','-r300');
exportgraphics(gcf,'transparentD.eps',...   % since R2020a
    'ContentType','vector',...
    'BackgroundColor','none')
close all