% Orginize feature data at a common "Cells" directory
% Input: features from the tasks directories (in the tracking folder)

% Assaf Zaritsky, Jan 2018
function [] = pcSingleCell_LEVER_LBP_SHAPE(params,dirs)

outdir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Cells/LEVER_LBP_SHAPE/';

if ~exist(outdir,'dir')
    error([outdir ' does not exist']);
end

cellTYXFname = [dirs.tracking 'cellIdTYX.mat'];
load(cellTYXFname);% cellTYX

lbpFname = [dirs.tracking filesep 'lbpDataLever.mat'];
load(lbpFname); % lbpData

shapeFname = [dirs.tracking filesep 'shapeDataLever.mat'];
load(shapeFname); % shapeData

if isempty(lbpData)
    return;
end

nCurCells = length(lbpData.fov);

LeverLbpShapeFname = [outdir filesep params.curFname sprintf('_s%02d',params.curTask) '_LEVER_LBP_SHAPE.mat'];

if exist(LeverLbpShapeFname,'dir')  && ~params.always
    return;
end

dataTask = cell(1,nCurCells);
for icell = 1 : nCurCells  
    dataTask{icell}.ts = cellTYX{icell}.ts;
    dataTask{icell}.xs = cellTYX{icell}.xs;
    dataTask{icell}.ys = cellTYX{icell}.ys;
    
    if isempty(lbpData.fov{icell})
        dataTask{icell} = [];
        continue;
    end
    
    dataTask{icell}.lbpFOV = lbpData.fov{icell}.lbp;
    dataTask{icell}.lbpFWD = lbpData.fwd{icell}.lbp;
    dataTask{icell}.lbpBCK = lbpData.bck{icell}.lbp;
    
    dataTask{icell}.dlbp10dFOV = lbpData.fov{icell}.dlbp10d;
    dataTask{icell}.dlbp10dFWD = lbpData.fwd{icell}.dlbp10d;
    dataTask{icell}.dlbp10dBCK = lbpData.bck{icell}.dlbp10d;  
    
    dataTask{icell}.dlbp1dFOV = lbpData.fov{icell}.dlbp1d;
    dataTask{icell}.dlbp1dFWD = lbpData.fwd{icell}.dlbp1d;
    dataTask{icell}.dlbp1dBCK = lbpData.bck{icell}.dlbp1d;
    
    dataTask{icell}.dlbp1dFOV_median = median(lbpData.fov{icell}.dlbp1d);
    dataTask{icell}.dlbp1dFWD_median = median(lbpData.fwd{icell}.dlbp1d);
    dataTask{icell}.dlbp1dBCK_median = median(lbpData.bck{icell}.dlbp1d);
    
    dataTask{icell}.corrFovFwd = lbpData.corr{icell}.fov_fwd;
    dataTask{icell}.corrFovBck = lbpData.corr{icell}.fov_bck;
    dataTask{icell}.corrFwdBck = lbpData.corr{icell}.fwd_bck;
    
    dataTask{icell}.shapeFeats = shapeData.shapeFeats{icell}.feats;    
end
save(LeverLbpShapeFname,'dataTask','nCurCells');
fprintf([params.curFname sprintf('_s%02d',params.curTask) '_LEVER_LBP_SHAPE.mat']);
end