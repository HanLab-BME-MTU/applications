% Script to obtain movie info from compTracks and calculate the density in 
% the last frame.
%
% Luciana de Oliveira, February 2017.

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/10/20161004/target';
saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170207/target';

%Define strings for directory hierarchy as needed
 rDDir = {'rD100'};%,'rD60','rD80','rD120','rD140','rD160'};
 aPDir = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
outDirNum =1:30;
lRDir = {'lR0p03'};%;'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'

%define intensity and space information

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: pay attention for the 'areaSideLen','timeStep',and 'sampleStep'
%values. For the low receptor densities 'areaSideLen=25' and for high densities
%it is 12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infoIntensitySpace = struct('intensityInfo', [1 0.3],'probDim',2,'areaSideLen',12);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
  
        %iterate through the different runs
    for outDirIndx = 1 : length(outDirNum)
                   
        
        %iterate through the different labeling ratios
    for lRDirIndx = 1 : length(lRDir)
  
                
                %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                
                %load compTracks
                tmp = load([currDir,filesep,'compTracks.mat']);
                tracksSim = tmp.compTracks;
                clear tmp
                
                %get movie info
                
   compTracksAggregState = genMovieInfoLastFrameFromTracksSparse(tracksSim);
   
                                             
                %save results
                
                saveDir = [saveRoot,filesep,rDDir{rDDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                
               mkdir(saveDir)
                
                save([saveDir,filesep,'compTracksAggregState'],'compTracksAggregState','-v7.3');
                
                clear compTracks
                
    end %for each labelRatio
            
    end %for each output
        
    end %for each aP
    
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
