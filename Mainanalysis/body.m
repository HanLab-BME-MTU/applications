
function [PROPERTIES,ero,labeled] = body(seg_img,coord,regmax,logihalo,plusminus,methodDeterm)       
% body  determins what areas of an image are occupied by cells and
%       calculates image properties
%
% SYNOPSIS       [PROPERTIES,ero,labeled] = body(seg_img,coord,regmax,logihalo,plusminus,methodDeterm)
%
% INPUT          seg_img : either original image or segmented image
%                          (depends on method)
%                coord : set of coordinates
%                regmax : binary image giving the areas of nuclei
%                logihalo : binary image giving the areas of halos
%                plusminus : distance a set of coordinates may have to an
%                            cell area and still belong to it
%                methodDeterm : 1 or 2. Says if clustering or image
%                               segmentation has been applied applied
%                               (changes what body actually does)
%
% OUTPUT         PROPERTIES : PROPERTIES(:,1)=coord(:,1);
%						 	  PROPERTIES(:,2)=coord(:,2);
%							  PROPERTIES(:,3)=belongsto(:);  (number of cluster - label)
%							  PROPERTIES(:,4)=numberOfOccurences(:);  (how many cells in the cluster
%							                                           this cell is in)
%							  PROPERTIES(:,5)=bodycount(:);  (area of the cluster with the number given in belongsto)
%							  PROPERTIES(:,6)=perimdivare(:);  (cluster)
%                ero : is the binary image of the areas occupied by cells                
%                labeled : bwlabeled ero
%
%
% DEPENDENCIES   body uses {nothing}
%                                  
%                body is used by { trackCells
%                                  testbutton}
%
% Colin Glass, Feb 04         






%f=figure;
%ero is a binary picture of the whole area occupied by cells
%belongsto 
bodyDiskSize=35;
background =[]; 
I2 =[]; 
I3 =[];
balevel =[];
bw =[]; 
ero=[];
labeled=[];
belongsto=[];
helpcoord=[];

helpcoord=round(coord);


%%%%%%%%%%%%%%%%%%%%%%%% Important%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%depending on which method of image analysis we use , there two different
%approaches

if methodDeterm==1
    ero =(seg_img==1 | seg_img==3);
    
elseif methodDeterm==2
    %first we kick out the background, 
	background = imopen(seg_img,strel('disk',bodyDiskSize));
	I2 = imsubtract(seg_img,background); 
	I3 = imadjust(I2, stretchlim(I2), [0 1]);
	balevel = graythresh(I3);
    
	ero = im2bw(I3,balevel)|regmax|logihalo; 
    
else
    error('body doesnt know which method to use (methodDeterm~= 1|2)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%fill out the holes
%ero = imfill(bw,'holes');
ero= imdilate(ero,strel('disk',8));
ero= imerode(ero,strel('disk',10));
ero= imdilate(ero,strel('disk',4));
ero= imerode(ero,strel('disk',3));
ero= imdilate(ero,strel('disk',2));
ero= imerode(ero,strel('disk',1));
ero= imdilate(ero,strel('disk',1));
ero= imdilate(ero,strel('disk',1));

%figure,imshow(ero),title('body');
labeled=bwlabel(ero);
belongsto=zeros(length(helpcoord),1);

%now we determin, to which group each set of coordinates belongs
for i=1:length(helpcoord);
        if labeled(helpcoord(i,2),helpcoord(i,1)) ~= 0
               belongsto(i)=labeled(helpcoord(i,2),helpcoord(i,1));
        else
                helper=[];
                where=[];
                img_h=[];
                img_w=[];
                x_1=[];
                x_2=[];
                y_1=[]; 
                y_2=[];
                
                [img_h,img_w]=size(seg_img);
                x_1= round(helpcoord(i,2) - plusminus);
				x_2= round(helpcoord(i,2) + plusminus);
				y_1= round(helpcoord(i,1) - plusminus);
				y_2= round(helpcoord(i,1) + plusminus);
                   
				if x_1 < 1
                      x_1 = 1;
				end 
				
				if x_2 > img_h
                      x_2 = img_h;    
				end 
				
				if y_1 <1
                      y_2 = y_2 - y_1 + 1;
                      y_1 = 1;
				end 
				
				if y_2 > img_w
                      y_2 = img_w;
				end 
                
                where=labeled(x_1:x_2,y_1:y_2);
                helper=find(where);
                if ~isempty(helper) 
                    
                        here=[];
                        noth=[];
                        uniIdx=[];
                        uniquehel=[];
                        numberOfOc=[];
                        
                        helper = sort(helper);
						[uniquehel, uniIdx] = unique(helper);
						%uniqueIdx returns the last occurence of the respective unique entry
						%having sorted m before, we can now count the number of occurences
						if size(uniquehel,1) > size(uniquehel,2);
                                uniIdx = [0;uniIdx];
						else
                                uniIdx = [0,uniIdx];
						end
                        numberOfOcc = [];
						numberOfOcc = diff(uniIdx); 
                        [noth,here]=max(numberOfOcc);
                        belongsto(i)=uniquehel(here);
                        
                        clear here;
                        clear noth;
                        clear uniIdx;
                        clear uniquehel;
                        clear numberOfOc;
                        
                        
                end
                clear helper;
                clear where;
                clear img_h
                clear img_w
                clear x_1
                clear x_2
                clear y_1
                clear y_2
                        
         end
        
 end

 onBackGround=find(belongsto==0);
 if ~isempty(onBackGround)
     coord(onBackGround,:)=[];
     belongsto(onBackGround)=[];
 end
 
%now we determin how many times a set of coordinates fall into the same
%labeled area,thus: how many nucloi per area
belo = sort(belongsto);
[uniqueEntries, uniqueIdx] = unique(belo);
%uniqueIdx returns the last occurence of the respective unique entry
%having sorted m before, we can now count the number of occurences
if size(uniqueEntries,1) > size(uniqueEntries,2);
        uniqueIdx = [0;uniqueIdx];
else
        uniqueIdx = [0,uniqueIdx];
end 
numberOfOccurences = diff(uniqueIdx); 

noSensibleInf=[];

bodycount=zeros(length(uniqueEntries),1);
for i=1:length(uniqueEntries);
        %area of a group devided by the amount of nucloi in the group. Average
        onebody = [];
        onebody = labeled==uniqueEntries(i);
        if length(find(onebody))>10
                areabod = length(find(onebody));
                perimdivare(i) = length(find(bwperim(onebody)))/areabod;
                bodycount(i) = areabod;
            else 
                
                noSensibleInf(end+1,1)=i;
                perimdivare(i) = 0;
                bodycount(i) = 0;
                
            end
        
end 



perimdivarea=sum(perimdivare)/length(perimdivare);
%one row of hell is equivalent to the properties of one group(area)
hell=zeros(length(uniqueEntries),4);
hell(:,1)=uniqueEntries(:);
hell(:,2)=numberOfOccurences(:);
hell(:,3)=bodycount(:);
hell(:,4)=perimdivare(:);

hell(noSensibleInf,:)=[];

    
%one row of PROPERTIES gives all information for one set of coordinates
PROPERTIES=zeros(length(coord),6);
PROPERTIES(:,1)=coord(:,1);
PROPERTIES(:,2)=coord(:,2);
PROPERTIES(:,3)=belongsto(:);

f=0;

noSensibleProp=[];
%fill the right information into the right spots of PROPERTIES, via belongsto(now PROPERTIES(:,3)) (PROPERTIES<->hell)
for i=1:length(coord);
        f=[];
        f= find (PROPERTIES(i,3)==hell(:,1));
        if ~isempty(f)
            PROPERTIES(i,4)=hell(f,2);
            PROPERTIES(i,5)=hell(f,3);
            PROPERTIES(i,6)=hell(f,4);
            
        else
            noSensibleProp(end+1,1)=i;
        end
        
end 

PROPERTIES(noSensibleProp,:)=[];
clear f;

%avarea is the average area of all cells
%avarea=sum(bodycount.*numberOfOccurences)/length(coord);
