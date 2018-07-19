sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20161114/rD4/aP0p5/lR0p2/reforTracks';
trackSimPath='/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112114/targetISruns/rD4/aP0p5';


%define reform parameters information
reformParam = struct('simTimeStep',0.01,'sampleTimeStep',0.1,'frameRange',[1 100],'pixelSize',0.09,'psfSigma',1);
%Define strings for directory hierarchy as needed
outDirNum = 1:10;
lRDir = {'lR0p2'};
%Define number of label ratio
numLabelRatio = length(lRDir);

tic
 for outDirIndx = 1 : length(outDirNum)
           fprintf('\nProcessing out = %',outDirNum(outDirIndx));
                   
             %For each label ratio, load the simulated compTracks        
            for lRDirIndx=1:numLabelRatio                          
                currOutDir = [trackSimPath,filesep,'out',int2str(outDirNum(outDirIndx)),filesep,lRDir{lRDirIndx}];
               tempCompTracks=load([currOutDir,filesep,'compTracks.mat']);
               compTracks=tempCompTracks.compTracks;                             
            
             %compare tracks 
             [ reformTrack] = reformTracks( compTracks,reformParam );
              save([sourceRoot,filesep,'reformTrack_',int2str(outDirNum(outDirIndx))],'reformTrack','-v7.3');
              elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',elapsedTime);
           end %for each labelRatio
 end
 

 
 