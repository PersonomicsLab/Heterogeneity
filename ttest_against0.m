load('imagingdata.mat', 'Zs', 'IDP_names') 

% run ttest for all 90 IDPs
for i = 1:size(Zs, 2)
    [~, p(i)] = ttest(Zs(:,i));
end
pFDR = mafdr(p);
keep = find(pFDR < 0.05);
Zsig = mean(Zs(:, keep));

% bar plot
x = floor(size(Zsig, 2));
figure
bar(Zsig)
ylabel('Mean Normative Deviation')
%%%% add imaging feature names
xticklabels(IDP_names(keep))
set(gca,'XLim',[0 (x+1)],'XTick',1:1:x)

% print(gcf, 'ttest_against0', '-dpng', '-r300')