%% false discovery rate correct
function [P,H,P_agg] = cellmtx_FDR(P,varargin)
H = P;
dims = size(P);     % assuming P is a cell matrix (i.e., that length(dims)=2)
thresh = 0.05;
if nargin > 1
    thresh = varargin{1};
end

P_agg = [];
for i=1:dims(1)
    for j=(i+1):dims(2)
        P_agg = [P_agg; P{i,j}(~isnan(P{i,j}))];
    end
end

[~,P_agg_fdr] = fdr(P_agg);
H_agg = P_agg_fdr < thresh;

for i=1:dims(1)
    for j=(i+1):dims(2)
        if P{i,j}
            idx = find(~isnan(P{i,j}));
            P{i,j}(idx) = P_agg_fdr(1:length(idx));
            H{i,j}(idx) = H_agg(1:length(idx));
            P_agg_fdr(1:length(idx)) = [];
            H_agg(1:length(idx)) = [];
        end
    end
end
end