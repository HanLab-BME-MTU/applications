%% get folder names
cond = {[pwd filesep 'Col-I'],[pwd filesep 'Col-V20'],[pwd filesep 'Col-V100']};
nameList = {'Col-I';'Col_V20';'Col_V100'};
figPath = pwd;
numCond = numel(cond);
%% variables
MS_integ = cell(numCond,1);
Pers_integ = cell(numCond,1);
DC_integ = cell(numCond,1);
%% run function
for iCond = 1:numCond
    [meanSpeedCell, persistenceCell, diffCoeffCell] = getMeanSpeedPersistenceDiffusionCoeffFromTraj(cond{iCond});
    % Make integrated cell
    MS_integ{iCond} = cell2mat(meanSpeedCell);
    Pers_integ{iCond} = cell2mat(persistenceCell);
    DC_integ{iCond} = cell2mat(diffCoeffCell);
end
h1 = figure;
h1.Position(3:4) = [200 150];
%% box plot: MS
boxPlotCellArray(MS_integ,nameList,1,0,1)
ylabel('Mean speed (um/min)')
title('Mean Speed')
hgexport(h1,strcat(figPath,'/MeanSpeed'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/MeanSpeed'),'-v7.3')
print(h1,strcat(figPath,'/MeanSpeed.tif'),'-dtiff')
%% box plot: Pers
boxPlotCellArray(Pers_integ,nameList,1,0,1)
ylabel('Persistence (1)')
title('Persistence')
hgexport(h1,strcat(figPath,'/Persistence'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/Persistence'),'-v7.3')
print(h1,strcat(figPath,'/Persistence.tif'),'-dtiff')
% 
% Persistence is a unitless measure. It is calculated as the average cosine
% of the angles between consecutive displacement vectors. This value ranges
% from -1 to 1:  
% 
% 1 indicates perfect persistence, where the movement is in a straight
% line. 
% 0 indicates completely random movement, with no directional persistence. 
% -1 indicates perfect anti-persistence, where the movement direction is
% completely reversed at each step. 

%% box plot: Diffusion Coefficient
boxPlotCellArray(DC_integ,nameList,1,0,1)
ylabel('Diff Coeff (um^2/min)')
title('Diffusion Coefficient')
hgexport(h1,strcat(figPath,'/Diffusion'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/Diffusion'),'-v7.3')
print(h1,strcat(figPath,'/Diffusion.tif'),'-dtiff')

