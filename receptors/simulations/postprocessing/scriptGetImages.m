sourceRoot = 'C:\Users\s169185\Documents\DATA';
saveDirImput='C:\Users\s169185\Documents\DATA';
%Define strings for directory hierarchy as needed
 rDDir = {'rD4'};
 aPDir = {'aP0p0'};
outDirNum = 1;
lRDir = {'lR0p4'};%'lR0p1','lR0p2',

%Define number of label ratio
numLabelRatio = length(lRDir);

% define parameters structures:
timeInfo = struct('simTimeStep',0.01,'sampleStep',0.1,'cutOffTime',80);
intensityInfo = struct('bgav',5000,'bgnoise',50,'scaleFactor',1000);
spaceInfo = struct('pixelSize',0.09,'psfSigma',1,'imsize',25);


fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
      
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        compTracksVec = cell(numLabelRatio,1);
        
        %Original output is not organized by label ratio since
        %receptorInfoLabeled for each label ratio is saved as a struct.
        %Iterate through the outputs, to pull out each receptorInfoLabeled
        %then the compTracks.  
        for outDirIndx = 1 : length(outDirNum)
            
                       
            %For each label ratio, the inner most directory, create the
            %directory and save compTracks.            
            for lRDirIndx=1:numLabelRatio
                %define save directory
                saveDir=[saveDirImput,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                filesep,lRDir{lRDirIndx},filesep,'images'];
                saveInfo = struct('saveVar',1,'saveDir',saveDir);
                
                fprintf('\n   Out = %d, lR = %s ',outDirIndx,lRDir{lRDirIndx});
                 %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                filesep,lRDir{lRDirIndx}];
                
             
                %Write compTracks
                tracksSim = load([currDir,filesep,'compTracks.mat']);
                
%               Get the images
                getImageStackFromTracksDirect_New(tracksSim,timeInfo,intensityInfo,spaceInfo,saveInfo)

                
                clear tracksSim
                
                fprintf('... done.');
                
            end %for each labelRatio
            
                       
        end %for each outDir
        
                                    
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
