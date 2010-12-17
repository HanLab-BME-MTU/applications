bin=0.1;

% edit scriptDetectGeneral.m;
scriptDetectGeneral;
scriptTrackGeneral;

[tracksInfoMat, tracksIndxMat] = convStruct2MatNoMS(tracksFinal);

%extract the number of columns in the Info Matrix, currently it is 8
[~,InfoCol]=size(tracksInfoMat);
[~,IndxCol]=size(tracksIndxMat);
nInfoCol=round(InfoCol/IndxCol);

%extract the x and y coordinates, respectively:
tracksMatxCord=tracksInfoMat(:,1:nInfoCol:end);
tracksMatyCord=tracksInfoMat(:,2:nInfoCol:end);

%calculate the x and y coordinates of the displacement:
displMatxCord=tracksMatxCord(:,2:end)-tracksMatxCord(:,1:end-1);
displMatyCord=tracksMatyCord(:,2:end)-tracksMatyCord(:,1:end-1);

%calculate the absolute value of the displacement:
displMatabs=sqrt(displMatxCord.^2+displMatyCord.^2);


%frameJtoJp1=1;  %if you enter j it will display the displ between j and j+1
%figure(1)
%hist(displMatxCord(:,frameJtoJp1), min(displMatxCord(:,frameJtoJp1)):bin:max(displMatxCord(:,frameJtoJp1)));
%figure(2)
%hist(displMatyCord(:,frameJtoJp1), min(displMatyCord(:,frameJtoJp1)):bin:max(displMatyCord(:,frameJtoJp1)));


%calculate the standard deviation:
[~,col]=size(displMatabs);
for i=1:col
    colVec=displMatabs(:,i);
    colVec(isnan(colVec))=[];
    meanDispl(i)=mean(colVec);
    stdDispl(i)=std(colVec);
    relError(i)=stdDispl(i)/meanDispl(i);
    
    figure(i)
    quiver(tracksMatxCord(:,i), tracksMatyCord(:,i), displMatxCord(:,i),displMatyCord(:,i))
    saveas(gcf,['spt-displ',int2str(i),'.fig']);
    
    figure(100+i)
    hist(displMatabs(:,i), min(displMatabs(:,i)):bin:max(displMatabs(:,i)));
    saveas(gcf,['spt-hist',int2str(i),'.fig']);
end

meanDispl
stdDispl
relError

save('spt-displ-stat.mat', 'meanDispl','stdDispl','relError','displMatxCord','displMatyCord','displMatabs')



% Now do the same for only the full tracks:
fullTracksDisplABS=displMatabs;
for i=1:col
    fullTracksDisplABS(isnan(fullTracksDisplABS(:,i)),:)=[];
end
%doesn't improve much:
meanFullTracksDispl=mean(fullTracksDisplABS)
stdFullTracksDispl=std(fullTracksDisplABS)
relFullTracksDispl=stdFullTracksDispl./meanFullTracksDispl

%calculate the cumulative sum:
cumDispl=cumsum(fullTracksDisplABS,2);
meanCumDispl=mean(cumDispl)
stdCumDispl=std(cumDispl)
relErrorCumDispl=stdCumDispl./meanCumDispl


break;
%Now extract the same information from the flow determined by the kymograph
%analysis. First browse into the folder where the flow is:

for i=1:9
    load(['flow0',int2str(i),'.mat']);
    
    %the flow between two pairs of frames can not be connected
    corrDisplx=flow(:,3)-flow(:,1);
    corrDisply=flow(:,4)-flow(:,2);
    corrDisplabs=sqrt(corrDisplx.^2+corrDisply.^2);
    
    figure(i)
    quiver(flow(:,1),flow(:,2),corrDisplx,corrDisply);
    saveas(gcf,['corr-displ',int2str(i),'.fig']);
    
    figure(100+i)
    hist(corrDisplabs, min(corrDisplabs):bin:max(corrDisplabs));
    saveas(gcf,['corr-hist',int2str(i),'.fig']);   
    clear flow
    
    corrDisplabs(isnan(corrDisplabs))=[]
    
    meanCorrDisplabs(i)=mean(corrDisplabs);
    stdCorrDisplabs(i)=std(corrDisplabs);
    relErrorCorrDisplabs(i)=stdCorrDisplabs(i)./meanCorrDisplabs(i);
end

meanCorrDisplabs
stdCorrDisplabs
relErrorCorrDisplabs





