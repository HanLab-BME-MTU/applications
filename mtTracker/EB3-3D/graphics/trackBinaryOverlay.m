function trackXY=trackBinaryOverlay(img,XLimit,YLimit,tracks,frameIdx,colorIndex,colormap)
%%
trackXY=img;
if(length(size(img))<3)
    trackXY=repmat(img,1,1,3);
end
    
tracksColors=colormap;

if(~isempty(tracks))
    sampling=100;
    for trIdx=1:length(tracks)
        
        t=tracks(trIdx);
        tIdx=find(t.f==(frameIdx));
        if(~isempty(tIdx))
 %%
            RGB=tracksColors(colorIndex(trIdx),:);
            X=t.x(1:tIdx); Y=t.y(1:tIdx); Z=t.z(1:tIdx);
                       
            X=X-XLimit(1);
            Y=Y-YLimit(1);          
            XRatio=size(img,2)/(XLimit(2)-XLimit(1));
            YRatio=size(img,1)/(YLimit(2)-YLimit(1));
            X=X*XRatio;
            Y=Y*YRatio;
            inIdx=(X>0)&(Y>0)&(X<=size(img,2))&(Y<=size(img,1));
            X=X(inIdx);
            Y=Y(inIdx);
 
            xSeg=max(1,round(linspaceNDim(Y(1:(end-1)),Y(2:(end)),sampling)));
            ySeg=max(1,round(linspaceNDim(X(1:(end-1)),X(2:(end)),sampling)));
            indx=sub2ind(size(trackXY),xSeg,ySeg,ones(size(xSeg)));
            trackXY(indx)=RGB(1);
            indx=sub2ind(size(trackXY),xSeg,ySeg,2*ones(size(xSeg)));
            trackXY(indx)=RGB(2);
            indx=sub2ind(size(trackXY),xSeg,ySeg,3*ones(size(xSeg)));
            trackXY(indx)=RGB(3);
        end
    end
    
end
