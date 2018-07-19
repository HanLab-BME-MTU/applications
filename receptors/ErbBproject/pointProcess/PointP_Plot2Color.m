
function img=PointP_Plot2Color(posA,posB,sf,A,FN)
%PointP_Plot2color takes two lists of corrdinates (x,y) with their corresponding
%uncertainies (dx,dy) and displays them in an image. Each point is represented as
%an ellipse centered at (x,y) with major and minor axes of (dx,dy).
%ImgSize is the size of the original image. sf is a multipicitive scaling
%factor; 10x is a good starting point. A is the gaussian amp, FN is the
%filename to be saved. sub is an optional parameter [xMin,xMax,yMin,yMax]
%that uses a sub section of the points
%
%Jeffrey L. Werbin
%Harvard Medical School
%
%Last Update: 12/10/2012
%
%


%posA is displayed as Green, posB is red

 %Sets image size to contain all points
 pos = [posA ; posB];
 pos = pos*sf;
 Low = [ min(min(pos(:,1))), min(min(pos(:,2)))];
 High = [max(max(pos(:,1))),max(max(pos(:,2)))];
 ImgSize = ceil([High(1)-Low(1),High(2)-Low(2)]);
 Size = ImgSize;
 
 pixNum=ceil(max(max(pos(:,3:4)))*5);
 if ~mod(pixNum,2)
    pixNum = pixNum+1;
 end;
 
 %sanitizes point process
 Low = Low/sf;
 High = High/sf;
 posA(:,1)=posA(:,1)-Low(1);
 posA(:,2)=posA(:,2)-Low(2);
 
 posB(:,1)=posB(:,1)-Low(1);
 posB(:,2)=posB(:,2)-Low(2);
 
 ImgSize2 = ImgSize+2*pixNum;
 
 IA = PointP_Plot2(posA,ImgSize2,sf,A,FN,pixNum); 
 IB = PointP_Plot2(posB,ImgSize2,sf,A,FN,pixNum); 
 
 %removes border
 IA = IA(pixNum:pixNum+Size(1),pixNum:pixNum+Size(2),:);
 IB = IB(pixNum:pixNum+Size(1),pixNum:pixNum+Size(2),:);
 

 %removes buffered edge
 IA = IA(pixnum:pixnum+ImgSize(1),pixnum:pixnum+ImgSize(2));
 IB = IB(pixnum:pixnum+ImgSize(1),pixnum:pixnum+ImgSize(2));
 
 img = IA;
 img(:,:,1)=IB(:,:,2);
 img(:,:,3)=img(:,:,3)+IB(:,:,3);

 %imwrite(img,FN,'TIFF');
 %hImage = imshow(img);
    
end

