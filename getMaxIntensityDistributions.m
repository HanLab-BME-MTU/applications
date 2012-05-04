function res = getMaxIntensityDistributions(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('lb', [3:10 11 16 21 41 61 81 101 141]);
ip.addOptional('ub', [3:10 15 20 40 60 80 100 140 200]);
ip.addParamValue('Tlft_f', 4); % lifetime threshold, in frames
ip.parse(varargin{:});
lb = ip.Results.lb;
ub = ip.Results.ub;

nc = numel(lb);
nd = numel(data);
% mCh = find(strcmp(data.source, data.channels));
mCh = 1;

res = struct([]);
for i = 1:nd
    %load([data(i).source 'Tracking' filesep 'trackAnalysis.mat']);
    load([data(i).source 'Tracking' filesep 'tracksPost.mat']);
    tracks = tracks([tracks.catIdx]==1);
    
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
        w = min(ip.Results.Tlft_f, round(lb(k)/data(i).framerate));
        res(i).maxA{k} = arrayfun(@(t) max(t.A(mCh,:)), itracks);
        res(i).maxA4{k} = arrayfun(@(t) max(t.A(mCh,1:w)), itracks);
        res(i).sum4{k} = arrayfun(@(t) sum(t.A(mCh,1:w)), itracks);
        res(i).int4{k} = arrayfun(@(t) cumsum(t.A(mCh,1:w)), itracks, 'UniformOutput', false);
        if ub(k)<=10
            res(i).int{k} = arrayfun(@(t) cumsum(t.A(mCh,:)), itracks, 'UniformOutput', false);
            res(i).sum{k} = arrayfun(@(t) sum(t.A(mCh,:)), itracks);
        end
        
        % lifetimes for given cohort
        res(i).lft{k} = lft(idx);
        
    end
    res(i).lft_all = lft;
    res(i).maxA_all = arrayfun(@(t) max(t.A(mCh,:)), tracks);
    res(i).maxA_all_pstd = arrayfun(@(t) max(t.A_pstd(mCh,:)), tracks);
end
