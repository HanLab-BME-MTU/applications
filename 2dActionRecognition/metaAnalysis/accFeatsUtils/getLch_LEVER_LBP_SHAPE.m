%% 
function [feats, cellID, TXY] = getLch_LEVER_LBP_SHAPE(dataTask)
n = length(dataTask);
feats = initFeats();
cellID = [];
TXY = {};
ncells = 0;
for icell = 1 : n
    
    if isempty(dataTask{icell})
        continue;
    end
    
    ncells = ncells + 1;
    
    cellID = [cellID, icell];
    TXY{ncells}.ts = dataTask{icell}.ts;
    TXY{ncells}.xs = dataTask{icell}.xs;
    TXY{ncells}.ys = dataTask{icell}.ys;
    
    feats.lbpFOV = [feats.lbpFOV,median(dataTask{icell}.lbpFOV,1)'];
    feats.lbpFOV = [feats.lbpFWD,median(dataTask{icell}.lbpFWD,1)'];
    feats.lbpFOV = [feats.lbpBCK,median(dataTask{icell}.lbpBCK,1)'];
    
    feats.dlbp10dFOV = [feats.dlbp10dFOV,median(dataTask{icell}.dlbp10dFOV,1)'];
    feats.dlbp10dFWD = [feats.dlbp10dFWD,median(dataTask{icell}.dlbp10dFWD,1)'];
    feats.dlbp10dBCK = [feats.dlbp10dBCK,median(dataTask{icell}.dlbp10dBCK,1)'];
    
    feats.dlbp1dFOV_median = [feats.dlbp1dFOV_median,dataTask{icell}.dlbp1dFOV_median];
    feats.dlbp1dFWD_median = [feats.dlbp1dFWD_median,dataTask{icell}.dlbp1dFWD_median];
    feats.dlbp1dBCK_median = [feats.dlbp1dBCK_median,dataTask{icell}.dlbp1dBCK_median];
    
    feats.corrFovFwd = [feats.corrFovFwd,dataTask{icell}.corrFovFwd];
    feats.corrFovBck = [feats.corrFovBck,dataTask{icell}.corrFovBck];
    feats.corrFwdBck = [feats.corrFwdBck,dataTask{icell}.corrFwdBck];
    
    feats.shapeFeats = [feats.shapeFeats, median(dataTask{icell}.shapeFeats,1)'];   
end
end
    


%%
function [] = initFeats()
feats.lbpFOV = [];
feats.lbpFWD = [];
feats.lbpBCK = [];

feats.dlbp10dFOV = [];
feats.dlbp10dFWD = [];
feats.dlbp10dBCK = [];

feats.dlbp1dFOV = [];
feats.dlbp1dFWD = [];
feats.dlbp1dBCK = [];

feats.dlbp1dFOV_median = [];
feats.dlbp1dFWD_median = [];
feats.dlbp1dBCK_median = [];

feats.corrFovFwd = [];
feats.corrFovBck = [];
feats.corrFwdBck = [];

feats.shapeFeats = [];
end
