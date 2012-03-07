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
        w = min(4, round(lb(k)/data(k).framerate));
        res(i).maxA{k} = arrayfun(@(t) max(t.A(1,:)), itracks);
        res(i).maxA4{k} = arrayfun(@(t) max(t.A(1,1:w)), itracks);
        res(i).sum4{k} = arrayfun(@(t) sum(t.A(1,1:w)), itracks);
        res(i).int4{k} = arrayfun(@(t) cumsum(t.A(1,1:w)), itracks, 'UniformOutput', false);
        if ub(k)<=10
            res(i).int{k} = arrayfun(@(t) cumsum(t.A(1,:)), itracks, 'UniformOutput', false);
            res(i).sum{k} = arrayfun(@(t) sum(t.A(1,:)), itracks);
        end
    end
    
    res(i).maxA_all = arrayfun(@(t) max(t.A(1,:)), tracks);
    res(i).maxA_all_pstd = arrayfun(@(t) max(t.A_pstd(1,:)), tracks);
end
