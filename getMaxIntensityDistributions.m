function res = getMaxIntensityDistributions(data, lb, ub)

nc = numel(lb);

nd = numel(data);
res = struct([]);
for i = 1:nd
    load([data(i).source 'Tracking' filesep 'trackAnalysis.mat']);
    
    
    for k = 1:nc
        lft = [tracks.lifetime_s];
        itracks = tracks(lb(k)<=lft & lft<=ub(k));
        % max intensities
        maxA = arrayfun(@(t) max(t.A(1,:)), itracks);
        res(i).maxA{k} = maxA;
    end
    
    res(i).maxA_all = arrayfun(@(t) max(t.A(1,:)), tracks);
    res(i).maxA_all_pstd = arrayfun(@(t) max(t.A_pstd(1,:)), tracks);
end
