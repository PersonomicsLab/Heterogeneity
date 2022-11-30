%% get the FAs of each cluster
for i = 1:6
%  meanZ1(:,i) = mean(Z{i}(kms_final{i}==1,:,1));
%  meanZ2(:,i) = mean(Z{i}(kms_final{i}==2,:,1));
end
FAname = [{'Corpus Collosum', 'ALIC', 'Corona Radiata', 'Cingulum Cingulate Gyrus', 'Cingulum Hippocampus', 'Superior Longitudinal Fasciculus', 'Uncinate Fasciculus'}];
GMVname = [{'Insular Cortex', 'Superior Temporal Gyrus', 'Middle Temporal Gyrus', 'Frontal Medial Cortex', 'Cingulate Gyrus', 'Frontal Orbital Cortex', 'Thalamus', 'Hippocampus', 'Amygdala', 'Ventral Striatum'}];
CTname = [{'Entorhinal', 'Fusiform', 'Inferiorparietal', 'Parahippocampal', 'Superior Frontal'}];
FCname = [{'DMN Amplitude', 'CEN Amplitude', 'Salience Network Amplitude', 'DMN/DMN Connectivity', 'DMN/CEN Connectivity', 'DMN/SN Connectivity', 'CEN/CEN Connectivity', 'CEN/SN Connectivity', 'SN/SN Connectivity'}];
position = {'A. ', 'B. ', 'C. ', 'D. ', 'E. ', 'F. '};
group = {'Anhedonia', 'Depressed Mood', 'Somatic Disturbance', 'Chronic', 'Late Onset', 'Acute Impairment'};
for i = 1:6
data = meanZ1;
FA(:,1) = [mean(data(1:3, i),1);  mean(data(4:5, i),1); mean(data(6:11, i),1); mean(data(12:13, i),1); mean(data(14:15, i),1); mean(data(16:17, i),1); mean(data(18:19, i),1)];
GMV(:,1) = [mean(data(21:22, i),1); mean(data(23:26, i),1); mean(data(27:32, i),1);  mean(data(33:34, i),1); mean(data(35:38, i),1);  mean(data(39:40, i),1);  mean(data(41:42, i),1); mean(data(43:44, i),1); mean(data(45:46, i),1); mean(data(47:48, i),1)];
CT(:,1) = [mean(data([49,55], i),1);  mean(data([50,56], i),1); mean(data([51,57], i),1); mean(data([52,58], i),1); mean(data([53, 59], i),1)];
FC(:,1) = [mean(data(61:63, i),1); mean(data(64:65, i),1); mean(data(66:67, i),1); mean(data([68,69,74], i),1); mean(data([70,71,75,76,79,80], i),1); mean(data([72,73,77,78,81,82], i),1); data(83, i); mean(data([84,85,86,87], i),1); data(88, i)];

data = meanZ2;
FA(:,2) = [mean(data(1:3, i),1);  mean(data(4:5, i),1); mean(data(6:11, i),1); mean(data(12:13, i),1); mean(data(14:15, i),1); mean(data(16:17, i),1); mean(data(18:19, i),1)];
GMV(:,2) = [mean(data(21:22, i),1); mean(data(23:26, i),1); mean(data(27:32, i),1);  mean(data(33:34, i),1); mean(data(35:38, i),1);  mean(data(39:40, i),1);  mean(data(41:42, i),1); mean(data(43:44, i),1); mean(data(45:46, i),1); mean(data(47:48, i),1)];
CT(:,2) = [mean(data([49,55], i),1);  mean(data([50,56], i),1); mean(data([51,57], i),1); mean(data([52,58], i),1); mean(data([53, 59], i),1)];
FC(:,2) = [mean(data(61:63, i),1); mean(data(64:65, i),1); mean(data(66:67, i),1); mean(data([68,69,74], i),1); mean(data([70,71,75,76,79,80], i),1); mean(data([72,73,77,78,81,82], i),1); data(83, i); mean(data([84,85,86,87], i),1); data(88, i)];

spiders(group{i}, position{i}, FA, FAname, GMV, GMVname, CT, CTname, FC, FCname)
end

%% function
function spiders(group, position, FA, FAname, GMV, GMVname, CT, CTname, FC, FCname)
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions/github_repo (1)')
tit = [position group];
mn = -1;
mx = .6;
figure('units','normalized','outerposition',[0 0 1 1]);
t = tiledlayout(4,1);
title(t, tit, 'Fontsize', 18, 'fontweight','bold')
nexttile 
spider_plot(FC(:,:)', 'AxesLabels', FCname, 'LabelFontSize',11, 'AxesPrecision', 2, 'AxesLimits', [mn, mn, mn, mn, mn, mn, mn, mn, mn; mx, mx, mx, mx, mx, mx, mx, mx, mx], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
title('Resting State fMRI')
set(gca,'FontSize',12)
nexttile
spider_plot(GMV(:,:)', 'AxesLabels', GMVname, 'LabelFontSize',11, 'AxesPrecision', 2, 'AxesLimits', [mn, mn, mn, mn, mn, mn, mn, mn, mn, mn; mx, mx, mx, mx, mx, mx, mx, mx, mx, mx], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none');
title('Gray Matter Volume')
set(gca,'FontSize',12)
nexttile
spider_plot(CT(:,:)', 'AxesLabels', CTname, 'LabelFontSize', 11, 'AxesPrecision', 2, 'AxesLimits', [mn, mn, mn, mn, mn; mx, mx, mx, mx, mx], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
title('Cortical Thickness')
set(gca,'FontSize',12)
nexttile
spider_plot(FA(:,:)','AxesLabels', FAname, 'LabelFontSize', 11, 'AxesPrecision', 2, 'AxesLimits', [mn, mn, mn, mn, mn, mn, mn; mx, mx, mx, mx, mx, mx, mx], 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none');
title('Fractional Anisotropy') 
set(gca,'FontSize',12)
lg  = legend({'Cluster 1','Cluster 2'}); 
lg.Layout.Tile = 'South';
end
