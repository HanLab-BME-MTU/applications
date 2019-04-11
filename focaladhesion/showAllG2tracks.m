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
    idG2 = idsClassifiedStruct.idGroup1;
    idG6 = idsClassifiedStruct.idGroup6;
end
%% Load tracksNA belonged to only G1
% adhAnalProc = faPackage.getProcess(7);
% tracksNAG2=adhAnalProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',indexG2);
initRiseProc = faPackage.getProcess(11);
[tracksNAG2,idxTracks]=initRiseProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',indexG2);
%% indexG2 filtering
if ~isempty(idxTracks)
    idG2 = false(size(idG2));
    idG2(indexG2(idxTracks))=true;
end

idG2 = getForceTransmittingG1(idG2,tracksNAG2);
indexG2  = find(idG2);
tracksNAG2=initRiseProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',indexG2);

% indexG2  = find(idG2)';
%% another filtering for redundant tracks
[~,idxFinalTracks,idOtherOverlappingTracks] = filterOverlappingTracks(tracksNAG2);
idG2_nonOverlapping = false(size(idG2));
idG2_nonOverlapping(indexG2(idxFinalTracks'))=true;
indexG2NO = find(idG2_nonOverlapping & idG2);
indexFinalTracks = find(idxFinalTracks);

tracksNAG2=initRiseProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',indexG2NO);

numTracksG2 = numel(tracksNAG2);
%% Load imgStack, forceStack and another stack if it exists.
[imgStack, tMap, imgStack2] = getAnyStacks(MD);
%% Go through each adhesion and save the figure
gPath = [faPackage.outputDirectory_ filesep 'tracksG2'];
if ~exist(gPath,'dir')
    mkdir(gPath)
end
iForceSlave = 1;
initOutFolder = fileparts(initRiseProc.outFilePaths_{2,iForceSlave});
initDataPath = [initOutFolder filesep 'data'];
nameTitle=['initialLag Class' num2str(2)];

for ii=1:numTracksG2
    curTrack = tracksNAG2(ii);
    IDtoInspect = indexG2NO(ii);
    OtherOverlappingIDs = indexG2(idOtherOverlappingTracks==indexFinalTracks(ii));
    additionalName = 'G2';
    [h, timeLagMasterAgainstForce,timeLagMasterAgainstMainSlave]=showSingleAdhesionTrackSummary(MD,curTrack,imgStack,tMap,imgStack2,IDtoInspect, gPath,additionalName);
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
                idG2(OtherOverlappingIDs)=true;
        end  
    else
        print(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.eps'),'-depsc2')
        print(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.tif'),'-dtiff')
        savefig(h,strcat(gPath,'/trackID',num2str(IDtoInspect),additionalName,'.fig'))
    end
    close(h)
    % Quantify timelag from curTrack2...
    initialRiseLagAgainstEachSlave{1}(ii,1) = timeLagMasterAgainstForce;
    initialRiseLagAgainstEachSlave{2}(ii,1) = timeLagMasterAgainstMainSlave;
end
% Retrieve the existing inital rise time delay
initRiseStruct = load([initDataPath filesep nameTitle '.mat'],'initialLagTogetherAdjusted','nameList2');   
initialLagTogetherAdjusted=initRiseStruct.initialLagTogetherAdjusted; %nameList2 contains the information on how the time lag was calculated.
nameList2 = initRiseStruct.nameList2;
% Update it and save it again
if numel(nameList2)<2
    initialLagTogetherAdjusted = initialRiseLagAgainstEachSlave(1);
else
    % identify the common slave
    if strcmp(nameList2{1}(end),'2') && strcmp(nameList2{2}(end),'2') %if the common slave is amp2 or ampTotal2,
        % rearrange accordingly
        initialLagTogetherAdjusted{1} = initialRiseLagAgainstEachSlave{1}-initialRiseLagAgainstEachSlave{2};
        initialLagTogetherAdjusted{2} = -initialRiseLagAgainstEachSlave{2};
    elseif strcmp(nameList2{1}(end),'g') && strcmp(nameList2{2}(end),'g') %if the common slave is forceMag,
        initialLagTogetherAdjusted{1} = -initialRiseLagAgainstEachSlave{1};
        initialLagTogetherAdjusted{2} = initialRiseLagAgainstEachSlave{2}-initialRiseLagAgainstEachSlave{1};
    else 
        initialLagTogetherAdjusted = initialRiseLagAgainstEachSlave;
    end
end
% Save it
nameTitleBackup=[nameTitle 'Original'];
movefile([initDataPath filesep nameTitle '.mat'],[initDataPath filesep nameTitleBackup '.mat'])
save([initDataPath filesep nameTitle '.mat'],'initialLagTogetherAdjusted','nameList2');   

if askUser
    idGroup1=idG2;
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