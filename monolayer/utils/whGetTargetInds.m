function [negCtrlInds,restInds,targetsInds,posCntrl] = whGetTargetInds(geneDayDiff,targetGenesStr)
nTargets = length(targetGenesStr);

nGeneDay = length(geneDayDiff);

negCtrlInds = [];
restInds = [];
targetsInds = cell(1,nTargets);
posCntrl.inds = [];
posCntrl.targetInds = cell(1,4);
% posCntrl.targetStrs = {'CDC42','RAC1','beta-PIX','RHOA'};
 
posCntrl.targetInds = cell(1,3);
posCntrl.targetStrs = {'CDC42','RAC1','beta-PIX'};

for i = 1 : nTargets
    targetsInds{i} = [];
end

for iGeneDay = 1 : nGeneDay  
    geneStr = geneDayDiff{iGeneDay}.geneStr;
    KD = geneDayDiff{iGeneDay}.KD;
    nameStr = nan;
    if KD == 0
       negCtrlInds = [negCtrlInds iGeneDay];
       continue;
    else if strcmp(geneStr,'CDC42') || strcmp(geneStr,'RAC1') || strcmp(geneStr,'beta-PIX')
            posCntrl.inds = [posCntrl.inds iGeneDay];
            if strcmp(geneStr,'CDC42')
                posCntrl.targetInds{1} = [posCntrl.targetInds{1} iGeneDay];
            else if strcmp(geneStr,'RAC1')
                    posCntrl.targetInds{2} = [posCntrl.targetInds{2} iGeneDay];
                else if strcmp(geneStr,'beta-PIX')
                        posCntrl.targetInds{3} = [posCntrl.targetInds{3} iGeneDay];
                        %                     else if strcmp(geneStr,'RHOA')
                        %                             posCntrl.targetInds{3} = [posCntrl.targetInds{4} iGeneDay];
                        %                         end
                    end
                end
            end
            continue;
        else
            check = false;
            for t = 1 : nTargets
                if strcmp(geneStr,targetGenesStr{t})
                    targetsInds{t} = [targetsInds{t} iGeneDay];
                    check = true;
                    break;
                end
            end            
            % > 50%
            if ~check
                restInds = [restInds iGeneDay];
            end
        end
    end    
end
end