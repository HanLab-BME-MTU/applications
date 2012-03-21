function res = getMaxIntensityDistributions(data, lb, ub)

nc = numel(lb);

nd = numel(data);
res = struct([]);
for i = 1:nd
    load([data(i).source 'Tracking' filesep 'trackAnalysis.mat']);
    
%     validGaps = arrayfun(@(t) max([t.gapStatus 4]), tracks)==4;
%     singleIdx = [tracks.nSeg]==1;
%     vis = [tracks.visibility];
%     idx_Ia = singleIdx & validGaps & vis==1;
%     tracks = tracks(idx_Ia);
    
    lft = [tracks.lifetime_s];
    for k = 1:nc
        idx = lb(k)<=lft & lft<=ub(k);
        itracks = tracks(idx);
        % max intensities        
        w = min(4, round(lb(k)/data(i).framerate));
        res(i).maxA{k} = arrayfun(@(t) max(t.A(1,:)), itracks);
        res(i).maxA4{k} = arrayfun(@(t) max(t.A(1,1:w)), itracks);
        res(i).sum4{k} = arrayfun(@(t) sum(t.A(1,1:w)), itracks);
        res(i).int4{k} = arrayfun(@(t) cumsum(t.A(1,1:w)), itracks, 'UniformOutput', false);
        if ub(k)<=10
            res(i).int{k} = arrayfun(@(t) cumsum(t.A(1,:)), itracks, 'UniformOutput', false);
            res(i).sum{k} = arrayfun(@(t) sum(t.A(1,:)), itracks);
        end
        
        % lifetimes for given cohort
        res(i).lft{k} = lft(idx);
        
    end
    res(i).lft_all = lft;
    res(i).maxA_all = arrayfun(@(t) max(t.A(1,:)), tracks);
    res(i).maxA_all_pstd = arrayfun(@(t) max(t.A_pstd(1,:)), tracks);
end
