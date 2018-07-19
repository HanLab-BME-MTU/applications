function [stats, key,flat] = collectStatsData()
    % find all MEF directories
    base = '/project/biophysics/jaqaman_lab/lamins/2015/20150602/';
    D = dir([base 'MEF*']);
    % load in saved statistics
    clear stats;
    for i=fliplr(1:length(D))
        stats(i) = load([base D(i).name filesep 'stats_20150622.mat']);
    end
    if(nargout > 1)
        key = [num2cell(1:32)',cellfun(@(ML) ML.movieListFileName_,{stats.ML},'Unif',false)'];
    end
    if(nargout > 2)
        fields = setdiff(fieldnames(stats),'ML');
        for fi = 1:length(fields)
            field = fields{fi};
            flat.(field) = [stats.(field)];
        end
        flat.totalChannels = length(flat.(field));
        flat.nChannels = cellfun('length',{stats.(field)});
        setStart = cumsum(flat.nChannels);
        binaryStart = false(1,length(flat.(field)));
        binaryStart(setStart) = true;
        flat.originalIndex = cumsum(binaryStart);
        flat.ML = {stats(flat.originalIndex).ML};
    end
end