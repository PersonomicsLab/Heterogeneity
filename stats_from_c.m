function [stats, ctrack, ptrack] = stats_from_c(stats_cell)

num_imag = length(stats_cell);
num_colors = max(stats_cell{1}(:,2));

stats = zeros(num_colors, num_colors, num_imag,  num_imag);

for i=1:num_imag
    for j=1:size(stats_cell{i},1)
        colors = stats_cell{i}(j,1:2);
        ctrack{i,j} = colors;
        pval = stats_cell{i}(j,6);
        ptrack{i,j} = pval;
        if pval < 0.05
            stats(colors(1), colors(2), i, i) = floor(-log10(pval));
        end
    end
end

end