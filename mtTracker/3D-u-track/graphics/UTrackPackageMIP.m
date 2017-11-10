function UTrackPackageMIP(UTrackPackage)
% Display a full MIP of the tracks processed by the package
% Possible TODO: 
%     - detections
%     - specifying a frame of reference
%     - specifying a colormap

MD=UTrackPackage.getOwner();
processProj=ExternalProcess(MD,'rawProj');

myColormap=uint8( ...
    [[0 255 255]; ... % blue "all" tracks
    [0 255 0]; ... % green mapped tracks
    [255 0 100]; ... % kinetochore tracks
    [255 0 200]; ... % kinetochore tracks
    ]);

tic
disp('MIP')
project1D(  MD, ...
            'name','fullMIPNoManifold','channelRender','grayRed', ...
            'processSingleProj',processProj, 'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99]);   
        
tracksFinal=UTrackPackage.getProcess(2).loadChannelOutput(1);
tracks=TracksHandle(tracksFinal);
toc;
tic;
disp('Overlay')
overlayProjTracksMovie(processProj,'tracks',tracks, ... 
            'colorIndx',ones(size(tracks)),'colormap',myColormap,'name','tracks');
toc;
