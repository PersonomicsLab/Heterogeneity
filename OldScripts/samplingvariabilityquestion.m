
    for n = 1:6
        index = randsample(1:length(Z{7}), length(Z{n}));
        temp{n} = Z{7}(index,:);
        avg{n} = mean(temp{n});
    end
%%
for n = 1:6
hetavg{n} = mean(mean(abs(temp{n})));
hetstd{n} = mean(std(temp{n}));
homavg{n} = mean(mean(abs(Z{n})));
homstd{n} = mean(std(Z{n}));
end

%%
for n = 1:6
    homo{n} = mean(abs(Z{n}));
    hetero{n} = mean(abs(temp{n}));
end

for n = 1:6
    [h(n),p(n)] = ttest(hetero{n},homo{n});
end