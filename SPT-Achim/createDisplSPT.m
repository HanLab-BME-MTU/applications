bin=0.1;
sdT=[];

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

%take only the full tracks (can be extended to closed tracks in the future):
fullTracksMatxCord=tracksMatxCord;
fullTracksMatyCord=tracksMatyCord;

[~,nCols]=size(fullTracksMatxCord);
for j=1:nCols
    badEntries=isnan(fullTracksMatxCord(:,j));
    fullTracksMatxCord(badEntries,:)=[];
    fullTracksMatyCord(badEntries,:)=[];
end


%calculate the x and y coordinates of the displacement:
displMatxCord=fullTracksMatxCord(:,2:end)-fullTracksMatxCord(:,1:end-1);
displMatyCord=fullTracksMatyCord(:,2:end)-fullTracksMatyCord(:,1:end-1);

%We assume here, that ONLY the FIRST image is the reference Frame. Then we
%have to take the cumulative sum:

cumDisplFieldMatx=cumsum(displMatxCord,2);
cumDisplFieldMaty=cumsum(displMatyCord,2);

for j=1:nCols-1
    displField(j).pos(:,1)=fullTracksMatxCord(:,1);
    displField(j).pos(:,2)=fullTracksMatyCord(:,2);
    
    %here it is assumed that the there exists as many transformations as displ    
    if ~isempty(sdT) && length(sdT)==length(displField)
        Tx=sdT(j,1);
        Ty=sdT(j,2);    
    else
        Tx=0;
        Ty=0;
    end

    %I think that they are just permutated:
    displField(j).vec(:,1)=cumDisplFieldMatx(:,j)+Ty;
    displField(j).vec(:,2)=cumDisplFieldMaty(:,j)+Tx;
end

    
save('displField.mat', 'displField');

for j=1:length(displField)
    figure(j)
    quiver(displField(j).pos(:,1),displField(j).pos(:,2),displField(j).vec(:,1),displField(j).vec(:,2))
    set(gca,'YDir','reverse','XTick',[],'YTick',[])
end




