function [] = overlayMovieAdhesionsAssemRate( MD )
%function [] = manualLabelingAfterInitialClassification( MD ) opens an
% interactive window with classified adhesion tracks so that the user can
% select several adhesions to change their assignment to a right class. 

% Show some instruction window to tell a user what to do

%% Load FA package
faPackage=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = faPackage.getProcess(8);
% p = classProc.funParams_;


%% Load tracksNA
adhAnalProc = faPackage.getProcess(7);
iChan = find(adhAnalProc.checkChannelOutput);
tracksNA=adhAnalProc.loadChannelOutput(iChan,'output','tracksNA');

%% Load imgStack, forceStack and anyother stack if it exists.
[imgStack, tMap, imgStack2] = getAnyStacks(MD);
% Launch pickAdhesion window with labeled adhesions with a right color and
% unlabed ones with white color. Get the right classes per newly selected
% adhesions
%% draw
% pickAdhesionTracksInteractive(tracksNA, imgStack,...
%     'movieData',MD,'tMap',tMap, 'imgMap2',imgStack2, 'drawingOnly',true,...
%     'Property','assemRate','PropRange','auto');
pickAdhesionTracksInteractive(tracksNA, imgStack,...
    'movieData',MD,'tMap',tMap, 'imgMap2',imgStack2, 'drawingOnly',true,...
    'Property','assemRate','PropRange',[0.00022 0.25]);
%% pick

pickAdhesionTracksInteractive(tracksNA, imgStack,...
    'movieData',MD,'tMap',tMap, 'imgMap2',imgStack2);
end

