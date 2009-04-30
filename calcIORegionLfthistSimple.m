function [segmentStatusVector,segmentEUdistVector] = calcIORegionLfthistSimple(lftInfo, maskPattern, mode);
% calculate lifetimes outside and inside segmented area
% 
% INPUT:    lftInfo 	= lifetime info 
%           maskPattern = BW mask of corresponding pattern 
%           
% OUTPUT:   segmentStatusVector = vector of length corresponding to the
%           number of trajectories found in lftInfo.Mat_lifetime, where the
%           status value is 1 if the majority of this trajectory's
%           positions are inside the pattern, and 0 otherwise
%
% Dinah Loerke, last modified 09/09/2008

[no,nf] = size(lftInfo.Mat_lifetime);

cmx = full(lftInfo.Mat_xcoord);
cmy = full(lftInfo.Mat_ycoord);

lmat = lftInfo.Mat_lifetime;
lvec = full(max(lmat,[],2));


%% masks defining IN and OUT areas
maskIN  = maskPattern;
maskOUT = ~maskPattern;


% euclidian distance map of the segmented image
EDimage = bwdist(logical(maskIN));
EDimage_out = bwdist(logical(maskOUT));


figure;
imshow(double(maskIN),[]);
hold on;

    
%% loop over all objects and determine each object's status as either
%% inside or outside the segmented region
for n=1:no
        
    % traj of this object
    upos = find( cmx(n,:)>0 & cmy(n,:)>0 );
    currx = full(cmx(n,upos));
    curry = full(cmy(n,upos));   
   
    if isempty(upos) 
        objectInStatus(n) = 0;
        objectEUdist(1,n) = nan;
        objectEUdist(2,n) = nan;
        continue
    end
        
     % object is defined as IN if the majority of its points are inside    
    % maskIN
    instatus = [];
    eudist = [];
    eudist_out = [];
    for k=1:length(upos)
        instatus(k) = maskIN(ceil(curry(k)),ceil(currx(k)));
        eudist(k) = EDimage(ceil(curry(k)),ceil(currx(k)));
        eudist_out(k) = EDimage_out(ceil(curry(k)),ceil(currx(k)));
    end
    
    if length(find(instatus))>(0.5*length(upos))
        objectInStatus(n) = 1;
    else
        objectInStatus(n) = 0;
    end
    objectEUdist(1,n) = nanmean(nonzeros(eudist));
    objectEUdist(2,n) = nanmean(nonzeros(eudist_out));
    
    
    if (objectInStatus(n) == 0)  
        plot(currx(1),curry(1),'b.');
    else
        plot(currx(1),curry(1),'r.');
    end
        
end % of for n-loop

segmentStatusVector = objectInStatus;
segmentEUdistVector = objectEUdist;

end % of function
    
    