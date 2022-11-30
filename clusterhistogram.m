figure;
hist(ka, 3)
xticks([2 3 4])
xlim([2 4])
xlabel('Optimal Cluster Number')
ylabel('Number of Bootstraps')
legend('Anhedonia', 'Low Mood', 'Somatic', 'Chronic', 'Late Onset', 'Acute Impairment')