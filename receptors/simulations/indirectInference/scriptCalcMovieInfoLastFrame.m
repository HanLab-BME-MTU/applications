% Script to obtain movie info from compTracks and calculate the density in 
% the last frame.
%
% Luciana de Oliveira, February 2017.

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2017/01/20170126/diffLabelRatio/probe';
saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/probe';

%Define strings for directory hierarchy as needed
rDDir  = {'rD20','rD40','rD60','rD80','rD100','rD120','rD140'};%,'rD20','rD40','rD60','rD80','rD100','rD120','rD140'
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%,'aP0p2','aP0p3','aP0p4','aP0p6','aP0p7','aP0p8'
outDirNum =1:30;
lRDir = {'lR0p26'};%
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
                
   detectionInfo = genMovieInfoLastFrameFromTracksSparse(tracksSim);
   
                                             
                %save results
                
                saveDir = [saveRoot,filesep,rDDir{rDDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                
%                mkdir(saveDir)
                
                save([saveDir,filesep,'detectionInfo'],'detectionInfo','-v7.3');
                
                clear compTracks
                
    end %for each labelRatio
            
    end %for each output
        
    end %for each aP
    
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
