function [BWmask1, BWmask2, loop1_Pix, loop2_Pix, interface_Pix, BWmask1_ext, loop1, loop2, interface]=intersecMaskPolygon(mask,polygon,plotResults)
%-------------------------------------------------------------
% Input:
%
% mask:     can either be a: m x 2 matrix containing the point coordinates 
%                            of the perimeter.
%                        or: m x n BWmatrix containing 0/1 with 1 on the
%                            perimeter.
% polygon:  n x 2 matrix containing the vertices of the polygon.
%-------------------------------------------------------------
% Output:
% BWmask1,2: 0/1 matrix, 1 indicates that this position is in the domain.
%            and BWmask2 >= BWmask1.
%
% !!! Important !!! The interface as well as the vertices (the start and 
% end points of the interface) belong completely to the larger domain:
% the BWmask2. This is an arbitrary choice but it ensures that the 
% no pix is assigned to more than one domain.


if nargin < 3 || isempty(plotResults)
    plotResults=0;
end

[rowsMask,colsMask]=size(mask);
if colsMask==2
    isCurve=1;
    % if the mask comes as a series of points and not as 0/1 matrix, then
    % create this matrix
    max_x=max(max(mask(:,1)),max(polygon(:,1)));
    max_y=max(max(mask(:,2)),max(polygon(:,2)));
    BWmask=zeros(max_y,max_x);
    indMask=sub2ind(size(BWmask),mask(:,2), mask(:,1));
    BWmask(indMask)=1;
    [rowsMask,colsMask]=size(mask);
else
    isCurve=1;
    % if the mask comes as a mask, make sure it has the right connectivity!
    maskStruct=bwboundaries(mask,4);
    if length(maskStruct)>1
        display('!!! Cell segmentation is not unique, redo cutOutForceField!!!');
        mask=[];
        for bdID=1:length(maskStruct)
            if length(maskStruct{bdID})>length(mask)
                mask=maskStruct{bdID};
            end
        end
    else
        mask=maskStruct{1};
    end
    max_x=max(max(mask(:,1)),max(polygon(:,1)));
    max_y=max(max(mask(:,2)),max(polygon(:,2)));
    BWmask=zeros(max_y,max_x);
    indMask=sub2ind(size(BWmask),mask(:,2), mask(:,1));
    BWmask(indMask)=1;    
    [rowsMask,colsMask]=size(mask);
end

%----------------
% begin old code
% pixelPath=[];
% for i = 1:size(polygon,1)-1
%     p1 = polygon(i  ,:);
%     p2 = polygon(i+1,:);
%     piece(i).pix = bresenham(p1,p2,4); % pList is a Mx2 matrix
%     if i==1
%         pixelPath=vertcat(pixelPath,piece(i).pix(1:end,:)); 
%     else
%         % To avoid double points:
%         pixelPath=vertcat(pixelPath,piece(i).pix(2:end,:));
%     end
% end
% pixelPath=removeDoublePoints(pixelPath);
% end old code
%----------------
% begin new code
[pixelPath,piece]=createPixelPath(polygon);
% begin new code
%----------------

indPixelPath = sub2ind(size(BWmask), pixelPath(:,2), pixelPath(:,1));    
indHits = find(BWmask(indPixelPath));

% Check if there are less then two points
if length(indHits)<2
    display('Sorry, something went wrong, please re-draw the interface')
    return;
else %There might be more than two points. Then take the ones with the 
     % largest distance:
    maxD=0;
    for i=1:length(indHits)-1
        cand1=pixelPath(indHits(i),:);
        cand2=pixelPath(indHits(i+1),:);
        dc1c2=sum((cand2-cand1).^2);        
        if dc1c2>maxD 
            xPoint1=cand1;
            xPoint2=cand2;
            maxD=dc1c2;
        end
    end    
end

% Determine the real interface start and end point beeing xPoint1 and
% xPoint2. To do this run through all interface pieces (or edges) and
% determine the edge containing the start and end point.
for i=1:length(piece)
    for j=1:length(piece(i).pix)
        if compPts(xPoint1,piece(i).pix(j,:))
            indxPt1=i;
        end
        if compPts(xPoint2,piece(i).pix(j,:))
            indxPt2=i;
        end        
    end
end
% reduce the polygon:
interface=polygon(indxPt1:indxPt2+1,:);
interface(1,:)=xPoint1;
interface(end,:)=xPoint2;


%find the position of these two points along the perimeter:
flippedL1=0;
flippedL2=0;
            
if isCurve==1
    loop1=[];
    loop2=[];
    fndxPT1=0;
    fndxPT2=0;
    for i=1:rowsMask
        if fndxPT1==0 && fndxPT2==0
            loop1=vertcat(loop1,mask(i,:));
            if compPts(mask(i,:),xPoint1) || compPts(mask(i,:),xPoint2)
                fndxPT1=1;
                % now starts the second loop:
                loop2=vertcat(loop2,mask(i,:));
                
                % fill in the interface for loop1
                if compPts(loop1(end,:),interface(1,:))
                    loop1=vertcat(loop1,interface(2:end,:));
                else                    
                    loop1=vertcat(loop1,flipud(interface(1:end-1,:)));
                    flippedL1=1;
                end
            end
        elseif fndxPT1==1 && fndxPT2==0
            loop2=vertcat(loop2,mask(i,:));
            if compPts(mask(i,:),xPoint1) || compPts(mask(i,:),xPoint2)
                fndxPT2=1;
                % now the first loop starts again:
                loop1=vertcat(loop1,mask(i,:));
  
                % fill in the interface for loop2
                if compPts(loop2(end,:),interface(1,:))
                    loop2=vertcat(loop2,interface(2:end,:));
                else
                    loop2=vertcat(loop2,flipud(interface(1:end-1,:)));
                    flippedL2=1;
                end
                % the second loop is now complete!
            end
        elseif fndxPT1==1 && fndxPT2==1
            loop1=vertcat(loop1,mask(i,:));
        end
    end    
end
% for numerical issues, the interface should be ordered in the same way, in
% both loops. The algorithm bresenham is not symmetric! Therefore 
% flip the loops back if they have been flipped:
if flippedL1==1
    loop1=flipud(loop1);
end
if flippedL2==1
    loop2=flipud(loop2);
end


% Also calculate the interface as integer position and determine the center 
% of the interface (for plotting):
interface_Pix.coord=[];
for i = 1:size(interface,1)-1
    p1 = interface(i  ,:);
    p2 = interface(i+1,:);
    interface_Pix.coord = vertcat(interface_Pix.coord,bresenham(p1,p2,4)); % pList is a Mx2 matrix
end
interface_Pix.coord=removeDoublePoints(interface_Pix.coord);
interface_Pix.center=interface_Pix.coord(ceil(length(interface_Pix.coord)/2),:);

% Calculate a 0/1 mask for each domain and determine which is the smallest.
% The smallest will be #1 the larger one will be #2.
loop1_Pix=[];
for i = 1:size(loop1,1)-1
    p1 = loop1(i  ,:);
    p2 = loop1(i+1,:);
    loop1_Pix = vertcat(loop1_Pix,bresenham(p1,p2,4)); % pList is a Mx2 matrix
end
loop1_Pix=removeDoublePoints(loop1_Pix);

loop2_Pix=[];
for i = 1:size(loop2,1)-1
    p1 = loop2(i  ,:);
    p2 = loop2(i+1,:);
    loop2_Pix = vertcat(loop2_Pix,bresenham(p1,p2,4)); % pList is a Mx2 matrix
end
loop2_Pix=removeDoublePoints(loop2_Pix);

BWmask1=zeros(max_y,max_x);
BWmask2=zeros(max_y,max_x);
indMask1=sub2ind(size(BWmask1),loop1_Pix(:,2), loop1_Pix(:,1));
indMask2=sub2ind(size(BWmask2),loop2_Pix(:,2), loop2_Pix(:,1));
BWmask1(indMask1)=1;
BWmask2(indMask2)=1;
% check if we have the right connectivity here:
BWmask1=imfill(BWmask1,'holes');
BWmask2=imfill(BWmask2,'holes');
% determine the number of pixel in each domain:
sizeBWmask1=sum(sum(BWmask1));
sizeBWmask2=sum(sum(BWmask2));
if sizeBWmask1>sizeBWmask2;
    storeLoop=loop1;
    loop1=loop2;
    loop2=storeLoop;
    
    storeLoop_Pix=loop1_Pix;
    loop1_Pix=loop2_Pix;
    loop2_Pix=storeLoop_Pix;    
    
    storeBWmask=BWmask1;
    BWmask1=BWmask2;
    BWmask2=storeBWmask;    
end
% This is the smaller mask with the interface, which we call the extended
% mask. For the BWmask2 these are the same!
BWmask1_ext=BWmask1;

% Now substract the interface from the smaller domain!
for i=1:length(interface_Pix.coord)
    BWmask1(interface_Pix.coord(i,2),interface_Pix.coord(i,1))=0;
end

% calculate the new boundaries:
loop1_Pix_old=loop1_Pix;
B1=bwboundaries(BWmask1,4,'noholes');
if length(B1)>1
    display('!!! Cell segmentation is not unique, redo cutOutForceField!!!');
    loop1_Pix=[];
    for bdID=1:length(B1)
        if length(B1{bdID})>length(loop1_Pix)
            loop1_Pix=fliplr(B1{bdID});
        end
    end
else
    loop1_Pix=fliplr(B1{1});
end

B2=bwboundaries(BWmask2,4,'noholes');
if length(B2)>1
    display('!!! Cell segmentation is not unique, redo cutOutForceField!!!');
    loop2_Pix=[];
    for bdID=1:length(B2)
        if length(B2{bdID})>length(loop2_Pix)
            loop2_Pix=fliplr(B2{bdID});
        end
    end
else
    loop2_Pix=fliplr(B2{1});
end



if plotResults==1
    figure(3)
    imagesc(abs(BWmask1-BWmask2*0.7))
    colormap('gray')
    hold on
    plot(loop1_Pix(:,1),loop1_Pix(:,2),'.g');
    plot(loop1_Pix_old(:,1),loop1_Pix_old(:,2),'xg');
    plot(loop2_Pix(:,1),loop2_Pix(:,2),'.b');   
    plot(interface_Pix.coord(:,1),interface_Pix.coord(:,2),'sb');
    plot(xPoint1(1),xPoint1(2),'or')
    plot(xPoint2(1),xPoint2(2),'or') 
    hold off
end

% Alternatively one could use the following approach, then the interface
% belongs partially to either cell.
% [PixPos_x PixPos_y]=meshgrid(1:max_x,1:max_y);
% BWmask1=inpolygon(PixPos_x,PixPos_y,loop1(:,1),loop1(:,2));
% BWmask2=inpolygon(PixPos_x,PixPos_y,loop2(:,1),loop2(:,2));
% 
% figure(1)
% imagesc(abs(BWmask1-0.7*BWmask2))
% hold on 
% plot(pixelPath(:,1),pixelPath(:,2),'.r')
% hold off

end
