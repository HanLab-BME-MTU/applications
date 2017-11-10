function functionMahalonobisDistanceMultipleParam(currFileProbe,currFileTarget,resultsDirectory,indexCell)

%Script to compare target and probe intermediate statistics data througth
%the calculation of Mahalanobis distance and pvalue. s and pValue matrix
% will be saved in same directory of the data in a file caled Results
%
% The outputs of this script are the sMatrix and pMatrix. They are 3D
% matrices where the rows are receptor densities, the colums are association probability
% and the third dimension is the label ratio.
%path where the cell array containing info about all the
%                   simulation and intermidiate statics parameters are saved.
%                   The cell array contains the folowing structure:
%                   column 1: file path;
%                   column 2: receptor density directory;
%                   column 3: dissociation rate directory;
%                   column 4: association probability directory;
%                   column 5: labeled fraction directory;
%                   column 6: receptor density values;
%                   column 7: association probability values;
%                   column 8: dissociation rate values;
%                   column 9: labeled fraction values;
%                   column 10: output values;
%Luciana de Oliveira, July 2016.
% Modifiade to deal with multiple parameters in May 2017, LRO.


% Load cell array with probe info

temp=load([currFileProbe,filesep,'cellInfoAllIntermidiateStatatistics.mat']);
cellInfoAllIntermidiateStatatisticsProbe=temp.cellInfoAllIntermidiateStatatistics;
sourceRootProbe=cellInfoAllIntermidiateStatatisticsProbe{1,1};
clear temp

%Define strings for directory hierarchy as needed
rDDir =cellInfoAllIntermidiateStatatisticsProbe{1,2};
dRDir =cellInfoAllIntermidiateStatatisticsProbe{1,3};
aPDir =cellInfoAllIntermidiateStatatisticsProbe{1,4};
lRDir =cellInfoAllIntermidiateStatatisticsProbe{1,5};


%load the cell file with intermediate statistics information for the target

temp=load([currFileTarget,filesep,'cellInfoAllIntermidiateStatatistics.mat']);
cellInfoAllIntermidiateStatatisticsTarget=temp.cellInfoAllIntermidiateStatatistics;
sourceRootTarget=cellInfoAllIntermidiateStatatisticsTarget{1,1};
clear temp

% for indexCell=1:size(cellInfoAllIntermidiateStatatisticsTarget,1)
rDtarget = cellInfoAllIntermidiateStatatisticsTarget{indexCell,2};
dRDirT=cellInfoAllIntermidiateStatatisticsTarget{indexCell,3};
aPtarget = cellInfoAllIntermidiateStatatisticsTarget{indexCell,4};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
lRtarget = cellInfoAllIntermidiateStatatisticsTarget{indexCell,5};


for rDTarIndx = 1 : length(rDtarget)
    
    for aPTarIndx = 1 : length(aPtarget)
        
        for dRDirIndx = 1 : length(dRDirT)
            
            for lRTarIndx = 1 : length(lRtarget)
                % load ratesAndDensity that has the paramMatrix
                currDir=[sourceRootTarget,filesep,rDtarget{rDTarIndx},filesep,...
                    dRDirT{dRDirIndx},filesep,aPtarget{aPTarIndx},filesep,lRtarget{lRTarIndx}];
                %                  currDir=[sourceRootTarget,filesep,rDtarget{rDTarIndx},filesep,...
                %                                  aPtarget{aPTarIndx},filesep,lRtarget{lRTarIndx}];
                tmp = load(fullfile(currDir,filesep,'ratesAndDensityComb_dt0p1_T10.mat'));
                
                % call paramMatrix and determine target cluster sizes
                
                paramMatrixTarget=tmp.paramMatrix;
                % clusterSizeTargetOnRate=size(tmp.rateOnPerClust,1);
                % clusterSizeTargetOffRate=size(tmp.rateOffPerClust,1);
                % clusterSizeTargetDensity=size(tmp.densityPerClust,1);
                
                %probes
                
                %,'rD6','rD8','rD10','rD12','rD14','rD16'};%,'rD4','rD6','rD8','rD10','rD12','rD14','rD16'};
                
                fprintf('\n===============================================================');
                
                
                %The top level directory is that of receptor density
                for lRDirIndx = 1 : length(lRDir)
                    
                    
                    tic
                    %Iterate through association probability values per density
                    for aPDirIndx = 1 : length(aPDir)
                        
                        fprintf('\nProcessing lR = %s, aP = %s ',lRDir{lRDirIndx},aPDir{aPDirIndx});
                        
                        %iterate through the different labeling ratios
                        for rDDirIndx = 1 : length(rDDir)
                            
                            for dRDirIndxP = 1 : length(dRDir)
                                currDir = [sourceRootProbe,filesep,rDDir{rDDirIndx},filesep,...
                                    dRDir{dRDirIndxP},filesep, aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
                                % % % %
                                %
                                %                   currDir = [sourceRootProbe,filesep,rDDir{rDDirIndx},filesep,...
                                %                            aPDir{aPDirIndx},filesep, dRDir{dRDirIndx},filesep,lRDir{lRDirIndx}];
                                
                                %   currDir = [sourceRootProbe,filesep,rDDir{rDDirIndx},filesep,...
                                %                 aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
                                
                                
                                % % % %
                                %name of current directory
                                %load ratesAndDensity that has the paramMatrix
                                tmp = load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                
                                %           call paramMatrix and determine probe cluster sizes
                                paramMatrixProbe=tmp.paramMatrix;
                                %                                 clusterSizeProbeOnRate=size(tmp.rateOnPerClust,1);
                                %                                 clusterSizeProbeOffRate=size(tmp.rateOffPerClust,1);
                                %                                 clusterSizeProbeDensity=size(tmp.densityPerClust,1);
                                
                                
                                %call function to compare target and probe data
                                [sMatrix(rDDirIndx,aPDirIndx,lRDirIndx),pMatrix(rDDirIndx,aPDirIndx,lRDirIndx)]=calculationMahalonobisDistanceandPvalue(paramMatrixTarget,paramMatrixProbe);
                                
                                %save file
                                
                                resultsDirBoot=[resultsDirectory,filesep,dRDirT{dRDirIndx},'Target',filesep,dRDir{dRDirIndxP},filesep,rDtarget{rDTarIndx},...
                                    aPtarget{aPTarIndx},lRtarget{lRTarIndx}];
                                
                                mkdir(resultsDirBoot)
                                save([resultsDirBoot,filesep,'sMatrix'],'sMatrix','-v7.3');
                                save([resultsDirBoot,filesep,'pMatrix'],'pMatrix','-v7.3');
                            end
                        end
                    end
                end
                
                %save s and p matrix
                %make a directory
                
            end
        end
    end
end
% end