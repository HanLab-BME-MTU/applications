function [] = showAllG1tracks( MD, askUser )
%function [] = showAllG1tracks( MD ) shows all G1 track summaries.
% input:    MD      MovieData
%           askUser     ask a user if the shown adhesion summary belongs
%                       to G1 class
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
idG1 = idsClassifiedStruct.idGroup1;
indexG1  = find(idG1)';
if askUser
    idG2 = idsClassifiedStruct.idGroup2;
    idG6 = idsClassifiedStruct.idGroup6;
    idG9 = idsClassifiedStruct.idGroup9;
end
%% Load tracksNA belonged to only G1
initRiseProc = faPackage.getProcess(11);
tracksNAG1=initRiseProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',indexG1);
%% indexG1 filtering
idG1 = getForceTransmittingG1(idG1,tracksNAG1);
indexG1  = find(idG1)';
tracksNAG1=initRiseProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',indexG1);

numTracksG1 = numel(tracksNAG1);
%% Load imgStack, forceStack and another stack if it exists.
[imgStack, tMap, imgStack2] = getAnyStacks(MD);
%% Go through each adhesion and save the figure
gPath = [faPackage.outputDirectory_ filesep 'tracksG1'];
if ~exist(gPath,'dir')
    mkdir(gPath)
end
for ii=1:numTracksG1
    curTrack = tracksNAG1(ii);
    IDtoInspect = indexG1(ii);
    additionalName = 'G1';
    h=showSingleAdhesionTrackSummary(MD,curTrack,imgStack,tMap,imgStack2,IDtoInspect);
    if askUser
        answer = listdlg('PromptString',['Does this belong to G1? ' 'Confirmation of each track' num2str(ii) '/' num2str(numTracksG1)],...
            'ListString',{'Yes(G1)','No(G6)','Others(G2)','Others(G9)'});
        % Handle response
        switch answer
            case 1
                if exist('gPath','var')
                    print(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.eps'),'-depsc2')
                    savefig(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.fig'))
                end
            case 2
                % take out the existing G2
                idG1(IDtoInspect)=false;
                % insert it to G6
                idG6(IDtoInspect)=true;
            case 3
                % take out the existing G2
                idG1(IDtoInspect)=false;
                % insert it to G6
                idG2(IDtoInspect)=true;
            case 4
                % take out the existing G2
                idG1(IDtoInspect)=false;
                % insert it to G6
                idG9(IDtoInspect)=true;
        end  
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
    idGroup9=idG9;
    
    save(classProc.outFilePaths_{4,iChan},'idGroup1','idGroup2','idGroup3','idGroup4',...
        'idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
end
disp('Done!')

end