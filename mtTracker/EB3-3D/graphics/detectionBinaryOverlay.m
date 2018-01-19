function detXY=detectionBinaryOverlay(img,XLimit,YLimit,detections,colorIndex,colormap,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('cumulative',false);
ip.addOptional('radius',2);
ip.parse(varargin{:});
p=ip.Results;
%%
detXY=img;
if(size(img,3)==1)
    detXY=repmat(img,1,1,3);
end
XRatio=size(img,2)/(XLimit(2)-XLimit(1));
YRatio=size(img,1)/(YLimit(2)-YLimit(1));
radius=ceil(p.radius*XRatio);
detColors=colormap;

if(~isempty(detections))
    for fIdx=1:length(detections)
        
        d=detections(fIdx);
        if(iscell(colorIndex))
            colIdx=colorIndex{fIdx};
        else
            colIdx=colorIndex;
        end
        if(~isempty(d)&&~(isempty(d.xCoord)))
            X=d.xCoord(:,1); Y=d.yCoord(:,1);% Z=t.z(1:tIdx);

            X=X-XLimit(1);
            Y=Y-YLimit(1);

            X=X*XRatio;
            Y=Y*YRatio;
            inIdx=(X>0)&(Y>0)&(X<=size(img,2))&(Y<=size(img,1));
            X=X(inIdx);
            Y=Y(inIdx);
            cIndex=colIdx(inIdx);
            drawingBoard=zeros(size(img,1),size(img,2));
            for dIdx=1:length(X)
                drawingBoard=MidpointCircle(drawingBoard,radius,Y(dIdx),X(dIdx),cIndex(dIdx));
            end
            uniqueCIdx=unique(cIndex);
            for ucIdx=1:length(uniqueCIdx)
                cIdx=uniqueCIdx(ucIdx);
                [I,J] = find(drawingBoard==cIdx);
                indx=sub2ind(size(detXY),I,J,1*ones(size(I)));
                detXY(indx)=colormap(cIdx,1);
                indx=sub2ind(size(detXY),I,J,2*ones(size(I)));
                detXY(indx)=colormap(cIdx,2);
                indx=sub2ind(size(detXY),I,J,3*ones(size(I)));
                detXY(indx)=colormap(cIdx,3);
            end
        end
    end

end

% Draw a circle in a matrix using the integer midpoint circle algorithm
% Does not miss or repeat pixels
% Created by : Peter Bone
% Created : 19th March 2007
function i = MidpointCircle(i, radius, xc, yc, value)


xc = int16(xc);
yc = int16(yc);
keeper=(xc<(size(i,1)-radius))&& ... 
       (yc<(size(i,2)-radius))&& ...
       (xc>radius)&& ...
       (yc>radius);
if(keeper)

x = int16(0);
y = int16(radius);
d = int16(1 - radius);

i(xc, yc+y) = value;
i(xc, yc-y) = value;
i(xc+y, yc) = value;
i(xc-y, yc) = value;

while ( x < y - 1 )
    x = x + 1;
    if ( d < 0 ) 
        d = d + x + x + 1;
    else 
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    i( x+xc,  y+yc) = value;
    i( y+xc,  x+yc) = value;
    i( y+xc, -x+yc) = value;
    i( x+xc, -y+yc) = value;
    i(-x+xc, -y+yc) = value;
    i(-y+xc, -x+yc) = value;
    i(-y+xc,  x+yc) = value;
    i(-x+xc,  y+yc) = value;
end
end
