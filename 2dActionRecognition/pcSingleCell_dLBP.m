% Measure: difference between LBP at consecutive time points
% Assaf Zaritsky, June 2016
function [] = pcSingleCell_dLBP(params,dirs)

outdir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Cells/dLBP/';

nScales = 4;
% scales = 1.0./2.^((1:nScales)-1);

cellTYXFname = [dirs.tracking 'cellIdTYX.mat'];
load(cellTYXFname);% cellTYX

lbpFname = [dirs.tracking 'lbpData.mat'];
load(lbpFname); % lbpData

nCurCells = length(lbpData.fov);

for iScale = 1 : nScales    
    
    curOutdir = [outdir filesep num2str(iScale) filesep ]; 
    if ~exist(curOutdir,'dir') 
        error([curOutdir ' does not exist']);
    end        
    
    dLbpFname = [curOutdir filesep params.curFname sprintf('_s%02d',params.curTask) '_dLBP.mat'];        
    
    if exist(dLbpFname,'dir')  && ~params.always
        continue;
    end        
    
    dLbpWell = cell(1,nCurCells);
    for icell = 1 : nCurCells
        ts = cellTYX{icell}.ts;
        xs = cellTYX{icell}.xs;
        ys = cellTYX{icell}.ys;
        lbp = lbpData.fov{icell}.pyramidLBP{iScale}.lbp;% 47 x 10
        dLbp = abs(lbp(1:end-1,:) - lbp(2:end,:));
        deltaLbp = sum(abs(dLbp),2);        
        
        dLbpWell{icell}.ts = ts;
        dLbpWell{icell}.xs = xs;
        dLbpWell{icell}.ys = ys;
        dLbpWell{icell}.lbp = lbp;
        dLbpWell{icell}.dLbp = dLbp;
        dLbpWell{icell}.deltaLbp = deltaLbp';        
        dLbpWell{icell}.deltaLbpMedian = median(deltaLbp); 
    end    
    save(dLbpFname,'dLbpWell','nCurCells');
    fprintf([params.curFname sprintf('_s%02d',params.curTask) '_dLBP.mat\n']);
end
end