function [] = showMostRepresentativeAssemblingAdhesion(MD)
%function [] = showMostRepresentativeAssemblingAdhesion(MD) searches
%through trajectories, pick one that shows bona fide new assembly and shows
%the best goodness of fit.
% input
%       MD          MovieData file that has run FAPackage
% output
%       figures will be generated and stored in
%       FocalAdhesionPackage/AssemblyFigure
% Sangyoon Han, March 4 2020

%% Let the user know whats going on
waitHan = msgbox({'You will see the cell window and classified adhesion tracks with color';
'label. After closing this window, you can select colored (classified)';
'adhesions to see their time-series and associated assembly rates'});

uiwait(waitHan);    
%% Load FA package
faPackage=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = faPackage.getProcess(8);
iChan = find(classProc.checkChannelOutput);

%% Load tracksNA
finalProc = faPackage.getProcess(11);
tracksNA=finalProc.loadChannelOutput(iChan,'output','tracksNA');
%% Load the idsClassifiedStruct
idsClassifiedStruct = finalProc.loadChannelOutput(iChan,'output','idClass');

%% Load imgStack, forceStack and anyother stack if it exists.
[imgStack, tMap, imgStack2] = getAnyStacks(MD);
%% Launch pickAdhesion window with labeled adhesions with a right color and
% unlabed ones with white color. Get the right classes per newly selected
% adhesions
[IDs, iGroups, iPrevGroups,tracksNA]=pickAdhesionTracksInteractive(tracksNA, imgStack,...
    'movieData',MD,'tMap',tMap, 'imgMap2',imgStack2, 'idSelected',idsClassifiedStruct);


