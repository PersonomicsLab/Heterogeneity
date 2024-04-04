%% aim 1 figure: Panel A-D
% load imaging data
load('imagingdata.mat', 'Z', 'Zs')
load('mainanalysis_results', 'sigidx', 'G')

% shorter names for the significant imaging features
fignames = {'Superior Longitudinal fasciculus FA (R)', ...
'Precentral Gyrus GMV (L)' , ...
'Precentral Gyrus GMV (R)', ...
'Insula GMV (L)' , ...
'Posterior Cingulate Gyrus GMV (L)' , ...
'Posterior Cingulate Gyrus GMV (R)' , ...
'Frontal Orbital Cortex GMV (L)' , ...
'Amygdala GMV (L)' , ...
'Amygdala GMV (R)' , ...
'Superiorfrontal CT (L)'};



%%%%%%Prep data and labels for bar plot
for i = 1:length(sigidx)
    [~, ~, stats] = anovan(Zs(:,sigidx(i)),G,'display','off');
   c{i} =  multcompare(stats);
end
close all
for i = 1:7
    meanZ(:,i) = mean(Z{i}(:,sigidx),1);
end
%%SEM
for i = 1:7
    for x = 1:10
        SEM(x,i) = std(Z{i}(:,sigidx(x)))/sqrt(length(Z{i}));
    end
end

%begin the figure
figure
tiledlayout(2,3)

%Panel A: GMV bar plot
nexttile([1 3])
meanZGMV = meanZ(2:9, :);
cGMV = c(2:9);
SEMGMV = SEM(2:9, :);
bar(meanZGMV);
ylabel('Normative Deviation Mean')
ylim([-.3 .3])
set(gca,'XLim',[0 9],'XTick',1:1:8)
xticklabels(fignames(2:9))
t = title('Gray Matter Volume (GMV)');
set(gca,'FontSize',14)
stat = stats_from_c(cGMV);
disp(['found total of ', num2str(length(find(stat))),' significant pairs'])
[stats_Y,stats_X1,stats_X2] = plot_stats(gca,stat); %https://www.mathworks.com/matlabcentral/fileexchange/92613-statistical-significance
t.FontSize = 16;
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(meanZGMV);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, meanZGMV(:,i), SEMGMV(:,i), 'k', 'linestyle', 'none');
end
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
vline([1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]);
%%%https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline?tab=reviews%2F2081773
legend({'Anhedonia','Low Mood','Somatic', 'Chronic','Late Onset', 'Acute Impairment', 'Heterogenous'},'Location','Northeast')

%Panel B: FA bar plot
nexttile
meanZFA = [meanZ(1, :); [0,0,0,0,0,0,0]];
empty = ones(21,6);
cFA = {c{1}, empty};
SEMFA = SEM(1, :);
bar(meanZFA);
ylabel('Normative Deviation Mean')
ylim([-.3 .3])
set(gca,'XLim',[0 2],'XTick',1:1:1)
xticklabels(fignames(1))
% t = title('Imaging Measures Significantly Different between Clinically Dissociated Groups');
t = title('Fractional Anisotropy (FA)');
set(gca,'FontSize',14)
stat = stats_from_c(cFA);
disp(['found total of ', num2str(length(find(stat))),' significant pairs'])
[stats_Y,stats_X1,stats_X2] = plot_stats(gca,stat); %https://www.mathworks.com/matlabcentral/fileexchange/92613-statistical-significance
t.FontSize = 16;
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(meanZFA);
ngroups = 1;
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, meanZFA(1,i), SEMFA(:,i), 'k', 'linestyle', 'none');
end
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off

%Panel C: CT bar plot
nexttile
meanZCT = [meanZ(10, :); [0,0,0,0,0,0,0]];
empty = ones(21,6);
cCT = {c{10}, empty};
SEMCT = SEM(10, :);
bar(meanZCT);
ylabel('Normative Deviation Mean')
ylim([-.3 .3])
set(gca,'XLim',[0 2],'XTick',1:1:1)
xticklabels(fignames(10))
t = title('Cortical Thickness (CT)');
set(gca,'FontSize',14)
stat = stats_from_c(cCT);
disp(['found total of ', num2str(length(find(stat))),' significant pairs'])
[stats_Y,stats_X1,stats_X2] = plot_stats(gca,stat); %https://www.mathworks.com/matlabcentral/fileexchange/92613-statistical-significance
t.FontSize = 16;
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(meanZCT);
ngroups = 1;
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, meanZCT(1,i), SEMCT(:,i), 'k', 'linestyle', 'none');
end
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off

% panel D: spider plot
nexttile
clear meanZ
for i = 1:7
 meanZ(:,i) = mean(Z{i},1);
end
mean_mods = [mean(meanZ(1:20,:),1); mean(meanZ(20:50,:),1); mean(meanZ(51:62,:),1); mean(meanZ(63:90,:),1)];
mod_names = [{'FA', 'GMV', 'CT', 'FC'}];
spider_plot(mean_mods(:,1:6)','AxesLabels', mod_names, 'LabelFontSize', 12, 'AxesPrecision', 3, 'AxesLabelsEdge', 'none', 'LabelFontSize', 12, 'AxesLimits', [-.13, -.13, -.13, -.13; .09, .09, .09, .09]); %, 'AxesDisplay', 'one', 'AxesLabelsEdge', 'none')
lg  = legend({'anhedonia','low mood','somatic', 'chronic','late onset', 'acute impairment'}); 

%print(gcf,'aim1_figure', '-dpng','-r300');