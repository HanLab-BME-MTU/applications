function [coor,regmax] = findNucloiTrack(newImg, lev, minsizenuc, maxsizenuc, methodDeterm)          
% findNucloiTrack detects dark areas and tries to fit cells into them
%
% SYNOPSIS       [coor,regmax]= findNucloiTrack(newImg, lev, minsizenuc, maxsizenuc, methodDeterm)
%
% INPUT          newImg : either original image or segmented image
%                          (depends on method)
%                lev : level used for minima detection
%                minsizenuc : minimal size for nuclei
%                maxsizenuc : maximal size for nuclei
%                methodDeterm : 1 or 2. Says if clustering or image
%                               segmentation has been applied (changes what findNucloiTrack actually does)
%
% OUTPUT         coor : found coordinates
%                regmax : binary image giving the areas of nuclei
%
% DEPENDENCIES   findNucloiTrack uses {nothing}
%                                  
%                findNucloiTrack is used by { ptTrackCells
%                                             testbutton}
%
% Colin Glass, Feb 04         


%newImg is the picture, or the segmented image (clustering), lev was determined from values taken from the picture
regmax = [];


%%%%%%%%%%%%%%%%%%%%%%%% Important%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%depending on which method of image analysis we use , there two different
%approaches

if methodDeterm==1
    regmax = newImg==1;
    
elseif methodDeterm==2
    regmax = imextendedmin(newImg,lev);
    
else
    error('findNucloiTrack doesnt know which method to use (methodDeterm~= 1|2)')
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
img_labels = [];

img_labels =  bwlabel(regmax);

huntdown = [];

huntdown = regionprops(img_labels,'Area','PixelList');

clear img_labels;


img_area = [];
destroy = [];

%find really small and big regions 
img_area = [huntdown.Area];
destroy = find( (img_area < round(minsizenuc/4*3)) | (img_area > maxsizenuc));
temp = huntdown(destroy);

clear img_area;
clear destroy;
clear huntdown;

%get all pixel within those regions
destcoord = cat(1,temp(:).PixelList);
clear temp;

%set those pixel to zero
for o = 1:size(destcoord,1)
	regmax(destcoord(o,2),destcoord(o,1)) = 0;
end

clear destcoord;



se = [];

%morphological operations
se = strel('disk',3);   
regmax = imclose(regmax,se);
SE = strel('disk',3); 
regmax = imdilate(regmax,SE);
clear se;

%label results
labeledone = bwlabel(regmax);
lastgroup = max(max(labeledone)) ;


%more morphological operations, but this time the regions are already labeled.
%So if one region gets seperated thus becoming two or more, it still has
%the same numbers.

Se = strel('disk',7);       % disk, radius 4
regmax = imerode(regmax,Se);
SE = strel('disk',6); 
regmax = imdilate(regmax,SE);



img_labels = [];
%now we label the image again. So we have one image labeled before the
%second morph oper and one labeled after
img_labels =  bwlabel(regmax);
counter = 0;
membership = [];

%for every group of labeledone (labeled before second morph oper) we allow
%one labeled group in img_labels (labeled after...). If there are more then 
%one region we selected the biggest.

for group = 1:lastgroup
	counter = counter+1;
	ranoutofvars = img_labels(find(labeledone==group));
	ranoutofvars = ranoutofvars(find(ranoutofvars));
	whichgroup = unique(ranoutofvars);
	if ~isempty(whichgroup)
		if length(whichgroup)>1.1  
			
			[uniqueEntries,numberOfOccurences] = countEntries(ranoutofvars);
			[useless,whichentry] = max(numberOfOccurences);
			membership(group,1) = uniqueEntries(whichentry);
		
		else 
			membership(counter,1) = whichgroup;
		
		end
	else
		counter = counter-1;
	end
end

%in membership we have at most one group (belonging to img_labels) for
%every each group belonging to labeledone
%regmax is overwritten now only including the groups (of img_labels) given in membership
if ~isempty(membership)
	membership = membership(find(membership));
	regmax = ismember(img_labels,membership);
	
	clear img_labels;
	img_labels = bwlabel(regmax);

end

Se = [];
clear Se;




realnucloi = [];
realnucloi = regionprops(img_labels,'Area','PixelList','MajorAxisLength','Centroid','Eccentricity');
%look up regionprops to add qualities of interest

clear img_labels;

%kick out really small and big areas
img_areare = [];
relevant = [];
tempor = [];

img_areare = [realnucloi.Area];
relevant = find(img_areare > minsizenuc);
tempor = realnucloi(relevant);

clear realnucloi;
clear img_areare;
clear relevant;

img_areare = [tempor.Area];
relevant = find(img_areare < maxsizenuc);
temporal = tempor(relevant);

clear tempor
clear img_areare;
clear relevant;


% % %meanecc = sum([temp2.Eccentricity])/size(temp2,1);


coor = round(cat(1,temporal.Centroid));



