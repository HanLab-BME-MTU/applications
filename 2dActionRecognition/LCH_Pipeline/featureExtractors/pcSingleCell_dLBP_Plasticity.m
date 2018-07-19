% Measure: plasticity
% Find sequences where first half of the sequence is different from second
% half (Wilcoxon test). Use 20 time-points.
% Assaf Zaritsky, May 2017
function [] = pcSingleCell_dLBP_Plasticity(params,dirs)

basedir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Cells/';
outdir = [basedir 'dLBP_Plasticity_t20_pval01/'];
dLbpDir = [basedir 'dLBP/'];

nScales = 4;

plasticityParams.timeInvl = 20;
plasticityParams.significance = 0.001;

for iScale = 1 : nScales    
    
    curOutdir = [outdir filesep num2str(iScale) filesep ]; 
    if ~exist(curOutdir,'dir') 
        error([curOutdir ' does not exist']);
    end        
    
    dLbpPlasticityFname = [curOutdir filesep params.curFname sprintf('_s%02d',params.curTask) '_dLBP_Plasticity.mat'];
    
    if exist(dLbpPlasticityFname,'dir')  && ~params.always
        continue;
    end  
    
    curDLbpDir = [dLbpDir filesep num2str(iScale) filesep ];
    dLBPFname = [curDLbpDir filesep params.curFname sprintf('_s%02d',params.curTask) '_dLBP.mat'];        
    load(dLBPFname); % dLbpWell{icell}.deltaLbp, nCurCells
              
    dLbpPlasticityWell = cell(1,nCurCells);
    for icell = 1 : nCurCells
        dLbpPlasticityWell{icell} = getDeltaLbpPlasticity(dLbpWell{icell}.deltaLbp,plasticityParams);
    end    
    save(dLbpPlasticityFname,'dLbpWell','nCurCells','dLbpPlasticityWell');    
end
end

%%
function dLbpPlasticity = getDeltaLbpPlasticity(deltaLbp,plasticityParams)

trajLength = length(deltaLbp);
plasticityTrajLength = trajLength - 2*plasticityParams.timeInvl + 1;

dLbpPlasticity.plasticityTrajLength = plasticityTrajLength;

dLbpPlasticity.wilcoxon = nan(1,plasticityTrajLength);

for i = 1 : plasticityTrajLength     
    curI = i + plasticityParams.timeInvl - 1;
    pval = ranksum(deltaLbp(curI-plasticityParams.timeInvl+1:curI),deltaLbp(curI+1:curI+plasticityParams.timeInvl));
    dLbpPlasticity.wilcoxon(i) = pval;
end

dLbpPlasticity.significant = dLbpPlasticity.wilcoxon < plasticityParams.significance;
dLbpPlasticity.nPlastic = countPlasticEvents(dLbpPlasticity.significant);

dLbpPlasticity.plasticity = dLbpPlasticity.nPlastic/dLbpPlasticity.plasticityTrajLength; % measure for normalized plasticity
assert(sum(isnan(dLbpPlasticity.wilcoxon)) == 0);
end

%%
function nPlastic = countPlasticEvents(significant)
nPlastic = 0;
i = 1;
while i < length(significant)
    if significant(i)
        nPlastic = nPlastic + 1;
        while (i <= length(significant))
            if significant(i)
                i = i + 1;
            else 
                break;
            end
        end
    end
    i = i + 1;
end
end