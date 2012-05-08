% Plots Point Process
%
%Jeffrey L. Werbin
%Harvard Medical School
%
%Last Update: 9/6/2011
%
%(posA,posB,sf,A,FN)
%

function PointP_Plot2Color(posA,posB,sf,A,FN)
%PointP_Plot2color takes two lists of corrdinates (x,y) with their corresponding
%uncertainies (dx,dy) and displays them in an image. Each point is represented as
%an ellipse centered at (x,y) with major and minor axes of (dx,dy).
%ImgSize is the size of the original image. sf is a multipicitive scaling
%factor; 10x is a good starting point. A is the gaussian amp, FN is the
%filename to be saved. sub is an optional parameter [xMin,xMax,yMin,yMax]
%that uses a sub section of the points

%posA is displayed as Green, posB is red

 %Sets image size to contain all points
 pos = [posA ; posB];
 pos = pos*sf;
 Low = [ min(min(pos(:,1))), min(min(pos(:,2)))];
 High = [max(max(pos(:,1))),max(max(pos(:,2)))];
 ImgSize = ceil([High(1)-Low(1),High(2)-Low(2)]);
 
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
 
 ImgSize = ImgSize+2*pixNum;
 
 IA = PointP_Plot2(posA,ImgSize,sf,A,FN,pixNum); 
 IB = PointP_Plot2(posB,ImgSize,sf,A,FN,pixNum); 

 img = IA;
 img(:,:,1)=IB(:,:,2);
 img(:,:,3)=img(:,:,3)+IB(:,:,3);

 imwrite(img,FN,'TIFF');
 hImage = imshow(img);
    
end

