function kinTracks=detectMTAlignment(MD,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('printAll',false, @islogical);
ip.addParameter('testKinIdx',[19 46 63 156],@isnumeric);
ip.addParameter('kinCapture',[]);
ip.addParameter('name','',@ischar);
ip.addParameter('bundlierPole',false);
ip.parse(MD,varargin{:});
p=ip.Results;
testKinIdx=p.testKinIdx;
printAll=p.printAll;

%%
if(isempty(p.kinCapture))
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'Kin' filesep 'catchingMT' filesep p.name];
    tmp=load([outputDirCatchingMT filesep 'catchingMT.mat'],'kinTracks');
    kinTracks=tmp.kinTracks;
else
    kinTracks=p.kinCapture;
end
% printAll=false;

%%  For each Kinetochore, indentify the connecting microtuble that bundle
% The pole ID use for the poleKin axis is defined by the location of the first point of each
% track
maxBundlingDistance=250;

for kIdx=1:length(kinTracks)
    progressText(kIdx/length(kinTracks),'Bundled MT.');

    kinTrack=kinTracks(kIdx);
    try
        kinTrack.addprop('bundledMatrix');       
        kinTrack.addprop('fiber');
    catch
    end;
    kinTrack.fiber=[];       
%     elevAziMin=nan(2,length(kinTrack.catchingMT));
%     elevAziMax=nan(2,length(kinTrack.catchingMT));
%     for mIdx=1:length(kinTrack.catchingMT)
%         mt=kinTrack.catchingMT(mIdx);
%         elevAziMin(:,mIdx)=[min(mt.elevation(mt.poleId(1),:)) min(mt.azimuth(mt.poleId(1),:)) ];
%         elevAziMax(:,mIdx)=[max(mt.elevation(mt.poleId(1),:)) max(mt.azimuth(mt.poleId(1),:)) ];
%     end
    KinXYStart=nan(2,length(kinTrack.catchingMT));
    KinXYEnd=nan(2,length(kinTrack.catchingMT));
    for mIdx=1:length(kinTrack.catchingMT)
        mtKinRef=kinTrack.catchingMTKinRef(mIdx);  
        KinXYStart(:,mIdx)=[(mtKinRef.x(1)) (mtKinRef.y(1)) ];
        KinXYEnd(:,mIdx)=[(mtKinRef.x(end)) (mtKinRef.y(end)) ];
    end
    
     % We know each trajectory finish in the kinetochore alignement 
     % Measuring the distance of starting point
     DStart = createSparseDistanceMatrix([KinXYStart'],[KinXYStart'], maxBundlingDistance);
     % Measuring the distance of the end point (tubular struture)
     DEnd = createSparseDistanceMatrix([KinXYEnd'],[KinXYEnd'], maxBundlingDistance);
     % Making sure the bundle is directed toward the kin
     % Naive version
     DStartEnd = createSparseDistanceMatrix([KinXYStart'],[KinXYEnd'], maxBundlingDistance);

     DStartEnd(logical(eye(size(DStartEnd))))=0;
     
     % For each microtubule collect the number of microtubule it is
     % connected to 
     
     bundledIdx=sum((DStart>0)&(DEnd>0)&(DStartEnd>0),1);
     kinTrack.fiber=full(bundledIdx);
     
     if(p.bundlierPole)
         mtPoleId=arrayfun(@(mt) mt.poleId(1), kinTrack.catchingMT);
         if(any((kinTrack.fiber>0)))
         [~,bundlierPoleId]=max(accumarray(mtPoleId(kinTrack.fiber>0),1));
         kinTrack.fiber=kinTrack.fiber(mtPoleId==bundlierPoleId);
         kinTrack.catchingMT=kinTrack.catchingMT(mtPoleId==bundlierPoleId);
         kinTrack.catchingMTKinRef=kinTrack.catchingMTKinRef(mtPoleId==bundlierPoleId);
         end
     end
%      bundledMatrix=((DStart>0)&(DEnd>0)&(DStartEnd>0));
%      kinTrack.bundledMatrix=bundledMatrix;
end

% First test, highlight bundle display the +TIP coordinate on a lateral view of the poleKin axis. 
outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'projections' filesep p.name filesep 'testBundleRadius' filesep]
system(['mkdir -p ' outputDirProj]);

if(printAll)    
    for kIdx=min(length(testKinIdx),testKinIdx)
        kinTrack=kinTracks(kIdx);
        
        [handles,~,fhandle]=setupFigure(1,2,'AxesWidth',8,'AxesHeight',4,'DisplayMode', 'print');
        hold(handles(1),'on');
        hold(handles(2),'on');
        
        for poleId=1:2
            scatter(handles(poleId),kinTrack.rho(poleId,:),zeros(size(kinTrack.rho(poleId,:))),'r');
            scatter(handles(poleId),0,0,'g');
        end
        
        for mIdx=1:length(kinTrack.catchingMT)
            mt=kinTrack.catchingMT(mIdx);
            mtKinRef=kinTrack.catchingMTKinRef(mIdx);
            
            %Project on the plan defined by the poleKin axis and the interpolar
            %axis.
            if(kinTrack.fiber(mIdx))
                plot(handles(mt.poleId(1)),mtKinRef.z,mtKinRef.x,'g-');
            else
                plot(handles(mt.poleId(1)),mtKinRef.z,mtKinRef.x,'b-');
            end
            ylim(handles(mt.poleId(1)),[-2000 2000])
            xlabel(handles(mt.poleId(1)),'Pole-Kinetochore axis (nm)')
            ylabel(handles(mt.poleId(1)),'Normal plane (nm)')
        end
             
        
        print([outputDirProj 'kin' num2str(kIdx,'%03d') '.png'],'-dpng');
        %%
        %     hold(handles(1),'off');
        %     hold(handles(2),'off');
        %close(fhandle);
    end
end

%%
outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles' filesep p.name];
system(['mkdir ' outputDirBundle]);
save([outputDirBundle filesep 'kin-MT-bundle.mat'],'kinTracks');

%% For test kinetochore, save and plot an Amira file with attached mt
outputDirAmira=[outputDirBundle filesep 'testKin' filesep 'Amira'];
system(['mkdir ' outputDirBundle  filesep 'testKin']);
%%

for kIdx=min(length(testKinIdx),testKinIdx)
    kinTrack=kinTracks(kIdx);
    trackSet=[kinTrack; kinTrack.catchingMT];
    trackType=[1; zeros(length(kinTrack.catchingMT),1)];
    bundleInfo=[0 kinTrack.fiber+1];
    amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',false,'edgeProp',{{'kinEB3',trackType},{'bundle',bundleInfo}})
end


if(p.printAll)
    %%

    %% For each kinetochore, plot an Amira file with attached mt
    outputDirAmira=[outputDirBundle filesep 'Amira_' filesep];
    parfor kIdx=1:length(kinTracks)
        kinTrack=kinTracks(kIdx);
        trackSet=[kinTrack; kinTrack.catchingMT];
        trackType=[1; zeros(length(kinTrack.catchingMT),1)];
        bundleInfo=[0 kinTrack.fiber+1];
        amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType},{'bundle',bundleInfo}})
    end
end


