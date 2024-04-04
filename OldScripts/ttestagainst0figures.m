 for n = 1:7
    [~, ttestP{n}, ~, stat] = ttest(Z{n}); %t-stat is basically just an inverst of p-value and the p-value for hetero group is of course larger, it's a much bigger group
tstat{n} = stat.tstat;
end

tstats = [tstat{1}, tstat{2}, tstat{3}, tstat{4}, tstat{5}, tstat{6}, tstat{7}];

%%
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data/allage')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Code')
 for n = 1:6
     [pthr,pcor,padj] = fdr(ttestP{n});
     sigkeepeach{n} = find(pcor < 0.05);
 end
 
 
 [~, pcortotal] = fdr(cell2mat(ttestP));
 ptotalshaped = reshape(pcortotal, 88, 7);
 for i = 1:7
 sigkeepall{i} = find(ptotalshaped(:,i) < 0.05);
 end
 
 for i = 1:7
 meanZ(:,i) = mean(Z{i},1);
 end

 for i = 1:7
    for x = 1:88
        SEM(x,i) = std(Z{i}(:,x))/sqrt(length(Z{i}));
    end
end
 
 %% visualize
titles = ["A. Anhedonia", "B. Low Mood", "C. Somatic Disturbance", "D. Chronic", "E. Late Onset", "F. Acute Impairment"];
names = ["Anhedonia", "Low Mood", "Somatic Disturbance", "Chronic", "Late Onset", "Acute Impairment"];
for i = 1:6
figure
bar(meanZ(:,i))
t = title(titles(i));
xticks([1:88])
xticklabels(IDPsselect89order)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'FontSize',8)
text(sigkeepall{i}, meanZ(sigkeepall{i},i), '*', 'FontSize',20);
t.FontSize = 16;
set(gcf, 'Position', [1 1 1312 976]);
view([90 90])
name = "ttestagainst0results" + names(i);
print(gcf,name,'-dpng','-r300');
end
 %%
 %%%run on big screen
 %%%then save it out with
 %print(gcf,'ttestagainst0results','-dpng','-r300');