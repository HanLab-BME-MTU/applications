function kMTComplex=createPairsComplex(kMTComplex,movieInfo) 


% Detectino mask parameters
theta=deg2rad(45);
D=10;
angleTest = @(x) x<=-cos(theta);

% figure
nComplex=numel(kMTComplex);
validFrames=find(arrayfun(@(x) ~isempty(x.xCoord),movieInfo))';
for iFrame=validFrames
    cometPos=[movieInfo(iFrame).xCoord(:,1) movieInfo(iFrame).yCoord(:,1)];
    cometAmp=movieInfo(iFrame).amp(:,1);
    
    kinPos1=arrayfun(@(x)x.kinetochores(1).pos(iFrame,:),kMTComplex,'Unif',false);
    kinPos1=vertcat(kinPos1{:});
    [index1,D1]=KDTreeBallQuery(cometPos,kinPos1,D);
    
    kinPos2=arrayfun(@(x)x.kinetochores(2).pos(iFrame,:),kMTComplex,'Unif',false);
    kinPos2=vertcat(kinPos2{:});
    [index2,D2]=KDTreeBallQuery(cometPos,kinPos2,D);

    %Debugging
%     figure; hold on;
%     plot(kinPos1(:,1),kinPos1(:,2),'^r');
%     plot(kinPos2(:,1),kinPos2(:,2),'^r');
%     plot(cometPos(:,1),cometPos(:,2),'og');
%     axis equal

    for iComplex=1:nComplex
        dP=kinPos2(iComplex,:)-kinPos1(iComplex,:);
        dP=dP/norm(dP);
        status1=false(numel(index1{iComplex}),1);
        for i=1:numel(index1{iComplex})            
            kP = cometPos(index1{iComplex}(i),:)-kinPos1(iComplex,:);
            status1(i)=angleTest(dot(kP,dP)/(norm(kP)));
        end
        
        if any(status1)
            mtIndex = index1{iComplex}(find(status1,1));
            kMTComplex(iComplex).microtubules(1).pos(iFrame,:)=cometPos(mtIndex,:);
            kMTComplex(iComplex).microtubules(1).amp(iFrame)=cometAmp(mtIndex);
        end
        
        status2=false(numel(index2{iComplex}),1);
        for i=1:numel(index2{iComplex})      

            kP = cometPos(index2{iComplex}(i),:)-kinPos2(iComplex,:);
            status2(i)=angleTest(-dot(kP,dP)/(norm(kP)));
        end
        
        if any(status2)
            mtIndex = index2{iComplex}(find(status2,1));
            kMTComplex(iComplex).microtubules(2).pos(iFrame,:)=cometPos(mtIndex,:);
            kMTComplex(iComplex).microtubules(2).amp(iFrame)=cometAmp(mtIndex);
        end
        
%         
%     for iTrack=find(validTracks)'
%         tracks1=tracks(pairsIndex(iTrack,1));
%         start1 = tracks1.seqOfEvents(1,1);
%         pairs(iFrame).pos1=tracks1.tracksCoordAmpCG(1+(iFrame-start1)*8:2+(iFrame-start1)*8);
% 
%         pairs(iFrame).pos2=tracks2.tracksCoordAmpCG(1+(iFrame-start2)*8:2+(iFrame-start2)*8);
%         
%         
% %         plot([pos1(1) pos2(1)],[pos1(2) pos2(2)],'o-r');
%         
%         kinPos(2*iTrack-1,:)=pos1;
%         kinPos(2*iTrack,:)=pos2;
%         
%         dP = pos2-pos1;
%         e = dP/norm(dP);
%         theta = atan2(dP(2),dP(1));
%         [xRange,yRange,nzIdx] = anisoGaussian2DSupport(pos1(1)-e(1)*kSigma*sigmaX/2,...
%             pos1(2)-e(2)*kSigma*sigmaX/2,sigmaX,sigmaY,theta,kSigma,MD.imSize_);
% %         plot([min(xRange) max(xRange) max(xRange) min(xRange) min(xRange)],...
% %             [min(yRange) min(yRange) max(yRange) max(yRange) min(yRange)],'-k');
%         [X,Y]=meshgrid(xRange,yRange);
%         ind=sub2ind(MD.imSize_,Y(nzIdx),X(nzIdx));
%         mask(ind)=true;
%         
%         [xRange,yRange,nzIdx] = anisoGaussian2DSupport(pos2(1)+e(1)*kSigma*sigmaX/2,...
%             pos2(2)+kSigma*e(2)*sigmaX/4,sigmaX,sigmaY,theta,kSigma,MD.imSize_);
% %         plot([min(xRange) max(xRange) max(xRange) min(xRange) min(xRange)],...
% %             [min(yRange) min(yRange) max(yRange) max(yRange) min(yRange)],'-k')
%         
%         [X,Y]=meshgrid(xRange,yRange);
%         ind=sub2ind(MD.imSize_,Y(nzIdx),X(nzIdx));
%         mask(ind)=true;
    end
end