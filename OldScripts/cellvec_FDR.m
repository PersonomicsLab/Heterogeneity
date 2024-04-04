function [p_cell,h,p_agg] = cellvec_FDR(p_cell, varargin)
h = p_cell;
n = length(p_cell);  % assumes p is 1D
p_agg = [];
thresh = 0.05;
if nargin > 1
    thresh = varargin{1};
end


% logical function for determining if a value is data or filler
cond_fun = @(x) (~isnan(x) & x~=0 & x~=1);
for i=1:n
    cellflag = false;
    data = p_cell{i};
    if iscell(data)
        cellflag = true;
        cdims = size(data); % assume length(cdims) <= 2
        if any(cdims == 1)
            [~,~,p_agg_i] = cellvec_FDR(data);
            p_agg = [p_agg; p_agg_i];
        else
            [~,~,p_agg_i] = cellmtx_FDR(data);
            p_agg = [p_agg; p_agg_i];
        end
    else
        p_agg = [p_agg; data(cond_fun(data))];
    end
end

[~,p_agg_fdr] = fdr(p_agg);
h_agg = p_agg_fdr < thresh;

if ~cellflag
    for i=1:n
        if ~(isempty(p_cell{i}))
            idx = find(cond_fun(p_cell{i}));
            p_cell{i}(idx) = p_agg_fdr(1:length(idx));
            h{i}(idx) = h_agg(1:length(idx));
            p_agg_fdr(1:length(idx)) = [];
            h_agg(1:length(idx)) = [];
        end
    end
end

end