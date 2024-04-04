% add paths
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Functions')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Data/allage')
addpath('/Users/kayla/Box Sync/ThesisLab_April2021/ThesisProjectCode/Code')
%load Z_TNM_case.csv
% Zs = table2array(readtable('Z_TNM_case.csv'));
load('finalPCAscores_kms_allage_Jun30.mat')
load('cellZ_new.mat')
kms = kms_final;
%load Z_TNM_case_eid.csv
% Zeids = table2array(readtable('Z_TNM_case_eid.csv'));

%%
%sort each cell in kms, then apply that sort to each cell in Z, then do
%Zlong = cell2mat(Z);, then do rho = corr(Zlong');
% for i = 1:6
%     [kmsort{i}, kmsidx{i}] = sort(kms{i});
% %     Zsort{i} = Z{i}(kmsidx{i}, :);
% end
% Zlong = [Zsort{1}; Zsort{2}; Zsort{3}; Zsort{4}; Zsort{5}; Zsort{6}];
% rho = corr(Zlong');

kmslong = [kmsort{1}; kmsort{2}; kmsort{3}; kmsort{4}; kmsort{5}; kmsort{6}];

%% visualize
%change the lines to the right numbers
figure; imagesc(rho); colorbar
hold on;
line([0,2436], [196, 196], 'Color', 'r');
line([0,2436], [313, 313], 'Color', 'r');
line([0,2436], [454, 454], 'Color', 'b');
line([0,2436], [592, 592], 'Color', 'b');
line([0,2436], [795, 795], 'Color', 'g');
line([0,2436], [928, 928], 'Color', 'g');
line([0,2436], [1105, 1105], 'Color', 'k');
line([0,2436], [1271, 1271], 'Color', 'k');
line([0,2436], [1586, 1586], 'Color', 'w');
line([0,2436], [1800, 1800], 'Color', 'w');
line([0,2436], [2156, 2156], 'Color', 'm');
line([0,2436], [2436, 2436], 'Color', 'm');

line([196, 196], [0,2436], 'Color', 'r');
line([313, 313], [0,2436], 'Color', 'r');
line([454, 454], [0,2436], 'Color', 'b');
line([592, 592], [0,2436], 'Color', 'b');
line([795, 795], [0,2436], 'Color', 'g');
line([928, 928], [0,2436], 'Color', 'g');
line([1105, 1105], [0,2436], 'Color', 'k');
line([1271, 1271], [0,2436], 'Color', 'k');
line([1586, 1586], [0,2436], 'Color', 'w');
line([1800, 1800], [0,2436], 'Color', 'w');
line([2156, 2156], [0,2436], 'Color', 'm');
line([2436, 2436], [0,2436], 'Color', 'm');

hold off
