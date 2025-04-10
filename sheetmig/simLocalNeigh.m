function []=simLocalNeigh

xLimit = [0 0];%[0  500];
yLimit = [0 0];%[0 2000];
numFrames = 240;
maxFrames = min(90,numFrames);
nCells = 1000;
Dx     = 1;
Dy     = 1;
xDrift = 0;
yDrift = 0;

[tracksMatxCord,tracksMatyCord]=simulateRandSheet(xLimit,yLimit,nCells,numFrames,Dx,Dy,xDrift,yDrift);

cell2edgeDist=zeros(size(tracksMatxCord));

%**************************************************************************
% copy paste part from cellMig_part_3_groupData:
%**************************************************************************

% Parameters for binning:
dFramesCorr=24;
maxDist=0;

% Parameters for plotting:
dframes=24;
dt=1/6;

% Don't store all the values, that would be crazy:
display('Measure local neighborhood:...')
for frame=1:dFramesCorr:maxFrames
    %         text=['Frame: ',num2str(frame,'%03.f'),' of ',num2str(numFrames)];
    %         progressText(frame/numFrames,text);
    
    % find the tracks which start at this frame
    xPos = tracksMatxCord(:,frame:numFrames); % Must be the same as tracksMatxCord because of cellDistClean (keep minTrackLength same)!!!
    yPos = tracksMatyCord(:,frame:numFrames); % Must be the same as tracksMatyCord because of cellDistClean (keep minTrackLength same)!!!
    c2eD =  cell2edgeDist(:,frame:numFrames);
    
    % Take only cells that are there at the start frame:
    cv1=~isnan(xPos(:,1));
    xPos=xPos(cv1,:);
    yPos=yPos(cv1,:);
    c2eD=c2eD(cv1,:);
    
    numCells=size(xPos,1);
    
    % This is to find the pairs of cell indices which correspond to
    % each row of pdist:
    % p is simply the ind in the vector:
    p=(1:(numCells-1)*numCells/2)';
    % vi is the index of "Spaltenumbruch" in the matrix:
    vi=1:(numCells-1);
    vi=vi*numCells-vi.*(vi+1)/2;
    %         try
    %             % find the corresponding indices with mex function:
    %             c2=comp2Vecs(p,vi');
    %         catch
    % The corresponding Matlab code which is five times slower than
    % the mex function:
    c2=zeros((numCells-1)*numCells/2,1);
    for i=1:length(p)
        c2(i)=sum(vi<p(i),2)+1;
    end
    %         end
    c1=p+c2.*(c2+1)/2-numCells*(c2-1);
    % The construction of c1 and c2 has been checked thoroughly!
    
    c2cDvec    =zeros(numCells*(numCells-1)/2,numFrames-frame+1);
    c2edCOM    =zeros(numCells*(numCells-1)/2,numFrames-frame+1);
    dirDisMean =zeros(numCells*(numCells-1)/2,numFrames-frame+1);
    trvDisMean =zeros(numCells*(numCells-1)/2,numFrames-frame+1);
    % travDisMean is the mean of the distances traveled by the pair of
    % cells.
    
    if numCells>1 % in case that there is at least 1 pair:
        for subframe=frame:numFrames
            % Calculate the cell to cell distance, the distance of the
            % center of mass to the wound edge and the direct distance and the traveled distance:
            c2cDvec(:,subframe-frame+1)= pdist([xPos(:,subframe-frame+1) yPos(:,subframe-frame+1)]);
            % c2edCOM(:,subframe-frame+1)= pdist( c2eD(:,subframe-frame+1),myMean);
            
            dirDist = sqrt((xPos(:,subframe-frame+1)-xPos(:,1)).^2+(yPos(:,subframe-frame+1)-yPos(:,1)).^2);
            
            trvDist = sum(sqrt((xPos(:,2:(subframe-frame+1))-xPos(:,1:(subframe-frame))).^2+(yPos(:,2:(subframe-frame+1))-yPos(:,1:(subframe-frame))).^2),2);
            
            %idVec=(1:numCells)';
            %c2 = pdist(idVec,retX)';
            %c1 = pdist(idVec,retY)';
            
            % [c1 c2]
            dirDisC1=dirDist(c1);
            dirDisC2=dirDist(c2);
            dirDisMean(:,subframe-frame+1)=mean([dirDisC1 dirDisC2],2);
            
            trvDisC1=trvDist(c1);
            trvDisC2=trvDist(c2);
            trvDisMean(:,subframe-frame+1)=mean([trvDisC1 trvDisC2],2);
            
            
            c2edC1=c2eD(c1,subframe-frame+1);
            c2edC2=c2eD(c2,subframe-frame+1);
            c2edCOM(:,subframe-frame+1)=mean([c2edC1 c2edC2],2);
        end
    else % there is not a single cell pair, might happen in single cell data:
        c2cDvec=NaN;
        c2edCOM=NaN;
        dirDisMean=NaN;
        trvDisMean=NaN;
    end
    % Take only the pairs which are there at the beginning and which
    % are close together:
    maxPairDist=1;
    cv1=~isnan(c2cDvec(:,1));
    cv2=c2cDvec(:,1)<maxPairDist;
    cvtot=cv1&cv2;
    
    c2cDvecOI    = c2cDvec(cvtot,:);
    c2edCOMOI    = c2edCOM(cvtot,:);
    dirDisMeanOI = dirDisMean(cvtot,:);
    trvDisMeanOI = trvDisMean(cvtot,:);
    
    % bin those according to the distance to the wound edge at the very
    % beginning:
    if maxDist==0
        bin=ones(size(c2edCOMOI(:,1)));
    else
        [~,bin]=histc(c2edCOMOI(:,1),0:bandWidth:maxDist+binPix);
    end
    
    if ~isempty(bin)
        for binID=1:max(bin)
            % find all cells that are in bin:
            ptsId=(bin==binID);
            % store the d(t) for each pair:
            neigh(binID,frame).c2cDis  =    c2cDvecOI(ptsId,:);
            neigh(binID,frame).c2cDD   =    c2cDvecOI(ptsId,:)-repmat(c2cDvecOI(ptsId,1),1,size(c2cDvecOI,2)); % change in relative distance over time
            neigh(binID,frame).dirDis  = dirDisMeanOI(ptsId,:);  % This is the mean value of the cell pair
            neigh(binID,frame).trvDis  = trvDisMeanOI(ptsId,:);  % This is the mean value of the cell pair
            
            % Also store the average distance to the wound edge:
            neigh(binID,frame).c2edMean = mean(c2edCOMOI(ptsId));
            neigh(binID,frame).c2edSTD  =  std(c2edCOMOI(ptsId));
        end
    else  % this happens for single cell migration if there is not a single pair
        neigh(1,frame).c2cDis  = [];
        neigh(1,frame).c2cDD   = [];
        neigh(1,frame).dirDis  = [];  % This is the mean value of the cell pair
        neigh(1,frame).trvDis  = [];  % This is the mean value of the cell pair
        neigh(1,frame).c2edMean = [];
        neigh(1,frame).c2edSTD  = [];
    end
    clear c2cDvec c2edCOM dirDisMean c2edMean;
end
groupData.neigh=neigh;
display('done');
% plotNeighResults(neigh,maxDist,dFramesCorr)
% close all



%**********************************************************************
% Perform measurement of local neigh:
%**********************************************************************
% put in by hand
marker=['ro','bs','m*','c+','gd','yo','ks'];

figure('Name','Local neighborhood');
clear M
label=1;
h=[];
%figure('Name',['Velocity correlation, for ',num2str(groupData.kPaClass(classID).yModu_kPa),' kPa-Class = ',groupData.kPaClass(classID).cc{1}],'NumberTitle','off')
% numFrames=groupData.glbMaxFrame(1);
numBins=1;
for binR=1:numBins
    numPlotFrames=length(1:dframes:maxFrames);
    k=1;
    for frame=1:dframes:maxFrames
        subplot(numPlotFrames,numBins,(k-1)*numBins+binR)
        
        % for classID=1:length(groupData.kPaClass)
        pixSize_um=0.64;
        
        % Calculate the mean separation speed:
        neigh=groupData.neigh(binR,frame);
        
        c2edClassMean    = nanmean(neigh.c2edMean);
        c2edSTDClassMean = nanmean(neigh.c2edSTD);
        
        % c2cDis c2cDD dirDis trvDis c2edMean c2edSTD
        
        % count paris in each bin:
        Ntot=sum(~isnan(neigh.c2cDD),1);
        c2cDDClassMean   =              nanmean(neigh.c2cDD,1);
        c2cDDClassSEM95  = facSEMtoSEM95*nanstd(neigh.c2cDD,[],1)./sqrt(Ntot);
        
        % Change in cell 2 cell distance, normalized by the average
        % distance traveled by the cell pair:
        normDDClassMean  =              nanmean(neigh.c2cDD./neigh.trvDis,1);
        normDDClassSEM95 = facSEMtoSEM95*nanstd(neigh.c2cDD./neigh.trvDis,[],1)./sqrt(Ntot);
        
        % c2cDDClassMean
        frameCount=length(c2cDDClassMean);
        % plot the results:
        %currh=errorbar(dt*((1:frameCount)-1),pixSize_um*c2cDDClassMean ,pixSize_um*c2cDDClassSEM95 ,['.',marker(2*(mod(classID,7)+1)-1)]);
        currh=errorbar(dt*((1:frameCount)-1),normDDClassMean,normDDClassSEM95);%,marker(2*(mod(classID,7)+1)-1:2*(mod(classID,7)+1)));
        h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
        % M{label}=[num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1},' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl),' ; <n>=',num2str(mean(Ntot),'%.1f')];
        label=label+1;
        hold on;
        plot(dt*((1:frameCount)-1),sqrt(2)./sqrt((1:frameCount)-1),'--k')
        plot(dt*((1:frameCount)-1),normDDClassMean-sqrt(2)./sqrt((1:frameCount)-1),'--r')
        % end
        title(['Cell to edge distance: ',num2str(pixSize_um*c2edClassMean,'%.1f'),'+-',num2str(pixSize_um*c2edSTDClassMean,'%.1f'),'[um]'])
        ylabel(['Change of cell to cell dist, dD(r,t=',num2str(dt*(frame-1)),'h)']);
        xlabel('dt [h]');
        ylim([ 0 0.4]);
        xlim([ 0 20]);
        %if frame ==1 && binR==numBins
        %    legend(h,M);
        %end
        % clear M
        label=1;
        h=[];
        hold off
        k=k+1;
    end
end
