function [] = showAllG2tracks( MD, askUser )
%function [] = showAllG1tracks( MD ) shows all G1 track summaries.
% input:    MD          MovieData
%           askUser     ask a user if the shown adhesion summary belongs
%                       to G2 class
% output:   if askUser is true, idClassifiedStruct is updated and
%           overwritten
%% Load FA package
faPackage=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = faPackage.getProcess(8);
% p = classProc.funParams_;
if nargin<2
    askUser = false;
end
%% Load the class ids
iChan = find(classProc.checkChannelOutput);
idsClassifiedStruct = load(classProc.outFilePaths_{4,iChan});
%% Get the G1 ids
idG2 = idsClassifiedStruct.idGroup2;
indexG2  = find(idG2)';
if askUser
    idG1 = idsClassifiedStruct.idGroup1;
    idG6 = idsClassifiedStruct.idGroup6;
end
%% Load tracksNA belonged to only G1
adhAnalProc = faPackage.getProcess(7);
tracksNAG2=adhAnalProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',indexG2);
%% indexG2 filtering
idG2_FT = getForceTransmittingG2(idG2,tracksNAG2,MD.timeInterval_);
% indexG2  = find(idG2)';
%% another filtering for redundant tracks
[~,idxFinalTracks,idOtherOverlappingTracks] = filterOverlappingTracks(tracksNAG2);
idG2_nonOverlapping = false(size(idG2));
idG2_nonOverlapping(indexG2(idxFinalTracks'))=true;
indexG2NO = find(idG2_nonOverlapping & idG2_FT);
indexFinalTracks = find(idxFinalTracks);

tracksNAG2=adhAnalProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',indexG2NO);

numTracksG2 = numel(tracksNAG2);
%% Load imgStack, forceStack and another stack if it exists.
[imgStack, tMap, imgStack2] = getAnyStacks(MD);
%% Go through each adhesion and save the figure
gPath = [faPackage.outputDirectory_ filesep 'tracksG2'];
if ~exist(gPath,'dir')
    mkdir(gPath)
end
for ii=1:numTracksG2
    curTrack = tracksNAG2(ii);
    IDtoInspect = indexG2NO(ii);
    OtherOverlappingIDs = indexG2(idOtherOverlappingTracks==indexFinalTracks(ii));
    additionalName = 'G2';
    h=showSingleAdhesionTrackSummary(MD,curTrack,imgStack,tMap,imgStack2,IDtoInspect, gPath,additionalName);
    if askUser
        answer = questdlg('Does this belong to G2?',['Confirmation of each track' num2str(ii) '/' num2str(numTracksG2)],...
                  'Yes(G2)','No(G6)','Others(G1)','No(G6)');
        % Handle response
        switch answer
            case 'Yes(G2)'
                if exist('gPath','var')
                    print(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.eps'),'-depsc2')
                    savefig(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.fig'))
                end
            case 'No(G6)'
                % take out the existing G2
                idG2(OtherOverlappingIDs)=false;
                % insert it to G6
                idG6(OtherOverlappingIDs)=true;
            case 'Others(G1)'
                % take out the existing G2
                idG2(OtherOverlappingIDs)=false;
                % insert it to G6
                idG1(OtherOverlappingIDs)=true;
        end  
    else
        print(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.eps'),'-depsc2')
        print(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.tif'),'-dtiff')
        savefig(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.fig'))
    end
    close(h)
end
if askUser
    idGroup1=idG1;
    idGroup2=idG2;
    idGroup3=idsClassifiedStruct.idGroup3;
    idGroup4=idsClassifiedStruct.idGroup4;
    idGroup5=idsClassifiedStruct.idGroup5;
    idGroup6=idG6;
    idGroup7=idsClassifiedStruct.idGroup7;
    idGroup8=idsClassifiedStruct.idGroup8;
    idGroup9=idsClassifiedStruct.idGroup9;
    
    save(classProc.outFilePaths_{4,iChan},'idGroup1','idGroup2','idGroup3','idGroup4',...
        'idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
end
disp('Done!')
end