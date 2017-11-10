function [ maskCrop,cropLimsX,cropLimsY,imgOut] = cropImageBasedOnMask( mask , pad,imgIn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin<3 
    imgIn = []; 
end 
[ny,nx] =size(mask); 
 idxPix =  find(mask);
    [yMask,xMask] = ind2sub([ny,nx],idxPix); 
    
  
    cropMinY = min(yMask)-pad;
    if cropMinY < 0
     
        cropMinY = min(yMask); 
    end 
        
    cropMinX = min(xMask)-pad;
    if cropMinX < 0
    
        cropMinX = min(xMask) ; 
    end 
        
    cropMaxX = max(xMask)+pad; 
    if cropMaxX > nx
     
        cropMaxX =max(xMask); 
    end 
    
    cropMaxY = max(yMask)+pad; 
    if cropMaxY > ny;
  
        cropMaxY = max(yMask);
    end 
   
   maskCrop = mask(cropMinY:cropMaxY,cropMinX:cropMaxX) ;
   if ~isempty(imgIn)
      imgOut =  imgIn(cropMinY:cropMaxY,cropMinX:cropMaxX); 
   end 
   cropLimsY = [cropMinY,cropMaxY]; 
   cropLimsX = [cropMinX,cropMaxX]; 


