%% Thhree channel experiment
repFolder='/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/2017_10_13/ChoK1_shRNA_WT_gTal_rVinc_frBead_010'; % With this I found the rep G2.

% repFolder='/storage/disk3/NA_recruitment/analysis/Sangyoon/NA_recruitment/2017_10_13/ChoK1_shRNA_WT_gTal_rVinc_frBead_003'; % For G1
% data loading
curMD = MovieData.load([repFolder filesep 'movieData.mat']);
iFAPack = curMD.getPackageIndex('FocalAdhesionPackage');
FAPack=curMD.packages_{iFAPack}; iTheOtherProc=9; iForceRead=10;
adhAnalProc = FAPack.processes_{7};
% tracksNAtal2all = adhAnalProc.loadChannelOutput(2,1);
tracksNAall=adhAnalProc.loadChannelOutput(2,'output','tracksNA');
idClasses=load(FAPack.processes_{8}.outFilePaths_{4,2});

% Now we have to combine this with readings from step 9 and 10
theOtherReadProc=FAPack.processes_{iTheOtherProc};
forceReadProc=FAPack.processes_{iForceRead};

if ~isempty(forceReadProc)
    forceReadObj = load(forceReadProc.outFilePaths_{1,2},'tracksForceMag'); % the later channel has the most information.
    tracksForceMag = forceReadObj.tracksForceMag;
    idxTracksObj = load(forceReadProc.outFilePaths_{2,2},'idxTracks');
    if ~isfield(idxTracksObj,'idxTracks')
        idxTracksObj = load(forceReadProc.outFilePaths_{6,2},'idxTracks');
    end
    idxTracks = idxTracksObj.idxTracks;
    tracksNAall2 = tracksNAall(idxTracks);
    if isfield(tracksForceMag,'forceMag')
        [tracksNAall2(:).forceMag] = tracksForceMag.forceMag;
    end
end

if ~isempty(forceReadProc)
    idGroup1 = idClasses.idGroup1(idxTracks);
    idGroup2 = idClasses.idGroup2(idxTracks);
    idGroup3 = idClasses.idGroup3(idxTracks);
    idGroup4 = idClasses.idGroup4(idxTracks);
    idGroup5 = idClasses.idGroup5(idxTracks);
    idGroup6 = idClasses.idGroup6(idxTracks);
    idGroup7 = idClasses.idGroup7(idxTracks);
    idGroup8 = idClasses.idGroup8(idxTracks);
    idGroup9 = idClasses.idGroup9(idxTracks);
else
    disp('Traction reading was not done. No further filtering...')
end
clear tracksNAall
tracksNAG2 = tracksNAall2(idGroup2);%tracksTalStruct.tracksG2;

tfmPack = curMD.packages_{1};
stracReadingProc = FAPack.processes_{10};
nFrames = curMD.nFrames_; tMap=[];
try
    T=load(tracReadingProc.outFilePaths_{3,1});
    T=T.T;
catch
    iBeadChan = 1; % might need to be updated based on asking TFMPackage..
    SDCProc_FA= FAPack.processes_{1};
    if ~isempty(SDCProc_FA)
        s = load(SDCProc_FA.outFilePaths_{3,iBeadChan},'T');    
        T_FA = s.T;
    else
        T_FA = zeros(nFrames,2);
    end

    SDCProc_TFM=tfmPack.processes_{1};
    %iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    if ~isempty(SDCProc_TFM)
        s = load(SDCProc_TFM.outFilePaths_{3,iBeadChan},'T');    
        T_TFM = s.T;
    else
        T_TFM = zeros(nFrames,2);
    end
    T = -T_TFM + T_FA;
end
for ii=nFrames:-1:1
    cur_tMap=tfmPack.processes_{4}.loadChannelOutput(ii,'output','tMap');
    cur_T = T(ii,:);
    cur_tMap2 = imtranslate(cur_tMap, cur_T(2:-1:1));
    tMap(:,:,ii) = cur_tMap2;
end
clear cur_tMap
sdcProc = FAPack.processes_{1};
imgStack = sdcProc.loadOutStack(2);
imgStack2 = sdcProc.loadOutStack(3);
%% See the iInit for these tracks
% disp(InitTimeLagMatureTal{iRepTal2})
% See the general trend of 8th track
% iNATal=10;
% showSingleAdhesionTrackSummary(MDtal,tracksNAtal2(iNATal),talImgStack,talForceStack,iNATal);
[naG2All,idsAll]=pickAdhesionTracksInteractive(tracksNAG2, imgStack, 'movieData',curMD,'tMap',tMap,'imgMap2',imgStack2);
% naG2All=[281   516   510   375   314   240   441   476   316];
%% Figure generation - G1 from G2 population
iNAG2= 136; % for cell3 newly analyzed. %568; %<-for cell3 %767;% <-for cell10 %465; %781; %1366; %781; % 3 is the one, and 7 is similar
showSingleAdhesionTrackSummary(curMD,tracksNAG2(iNAG2),imgStack,tMap,imgStack2,iNAG2); %naG2All{iRepTal2}(iNATal));
% Save the figure
%% Figure generation - G1
iNAG2= 179; % for cell3 newly analyzed. %568; %<-for cell3 %767;% <-for cell10 %465; %781; %1366; %781; % 3 is the one, and 7 is similar
showSingleAdhesionTrackSummary(curMD,tracksNAG2(iNAG2),imgStack,tMap,imgStack2,iNAG2); %naG2All{iRepTal2}(iNATal));
% Save the figure

%% G1
tracksNAG1 = tracksNAall2(idGroup1);%tracksTalStruct.tracksG2;
%% See the iInit for these tracks
[naG1All,idsAllG1]=pickAdhesionTracksInteractive(tracksNAG1, imgStack, 'movieData',curMD,'tMap',tMap,'imgMap2',imgStack2);
%% Figure generation
iNAG1= 2768; %2924; %465; %781; %1366; %781; % 3 is the one, and 7 is similar
showSingleAdhesionTrackSummary(curMD,tracksNAG1(iNAG1),imgStack,tMap,imgStack2,iNAG1); %naG2All{iRepTal2}(iNATal));
% I couldn't find the right G1 rep in 2017_10_13/ChoK1_shRNA_WT_gTal_rVinc_frBead_010. I'll try with other cell
