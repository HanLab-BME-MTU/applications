function [coor,regmax]= findnucloitrack(info,lev,minsizenuc,maxsizenuc)



%info is the picture, lev was determined from values taken from the picture
regmax=[];

regmax=imextendedmin(info,lev);
%figure, imshow(regmax);


   img_labels=[];
   
img_labels= bwlabel(regmax);

   huntdown=[];
   
huntdown=regionprops(img_labels,'Area','PixelList');%%%%%,'MajorAxisLength','Centroid','Eccentricity');
%look up regionprops to add qualities of interest

   clear img_labels;


   img_area=[];
   destroy=[];
   
%kick out really small and big areas
img_area=[huntdown.Area];
destroy=find(or(img_area < round(minsizenuc/4*3),img_area > maxsizenuc));
temp=huntdown(destroy);

   clear img_area;
   clear destroy;
   clear huntdown;

destcoord=cat(1,temp(:).PixelList);
   clear temp;

for o=1:size(destcoord,1)
       regmax(destcoord(o,2),destcoord(o,1))=0;
end

   clear destcoord;

   se=[];
   
se = strel('disk',3);   
regmax=imclose(regmax,se);
SE = strel('disk',3); 
regmax= imdilate(regmax,SE);
   clear se;

labeledone=bwlabel(regmax);

lastgroup=max(max(labeledone)) ;


   
Se = strel('disk',7);       % disk, radius 4
regmax = imerode(regmax,Se);
SE = strel('disk',6); 
regmax= imdilate(regmax,SE);



  img_labels=[];
   
img_labels= bwlabel(regmax);
counter=0;
for group=1:lastgroup
     counter =counter+1;
     ranoutofvars=img_labels(find(labeledone==group));
     ranoutofvars=ranoutofvars(find(ranoutofvars));
     whichgroup=unique(ranoutofvars);
     if ~isempty(whichgroup)
          if length(whichgroup)>1.1  
             
              [uniqueEntries,numberOfOccurences] = countEntries(ranoutofvars);
              [useless,whichentry]=max(numberOfOccurences);
              membership(group,1)=uniqueEntries(whichentry);
             
          else 
             membership(counter,1)= whichgroup;
	
          end
    else
          counter=counter-1;
    end
end

membership=membership(find(membership));
regmax=ismember(img_labels,membership);

 clear img_labels;
img_labels= bwlabel(regmax);

   Se=[];

   clear Se;
  



   realnucloi=[];
realnucloi=regionprops(img_labels,'Area','PixelList','MajorAxisLength','Centroid','Eccentricity');
%look up regionprops to add qualities of interest

 clear img_labels;
 
%kick out really small and big areas
   img_areare=[];
   relevant=[];
   tempor=[];
   
img_areare=[realnucloi.Area];
relevant=find(img_areare > minsizenuc);
tempor=realnucloi(relevant);

   clear realnucloi;
   clear img_areare;
   clear relevant;
   
img_areare=[tempor.Area];
relevant=find(img_areare < maxsizenuc);
temporal=tempor(relevant);

   clear tempor
   clear img_areare;
   clear relevant;



% % 
% % clear img_area;
% % img_area=[temp.Area];
% % relevant2=find(img_area < 3000);
% % temp2=temp(relevant2);
% % %meanecc=sum([temp2.Eccentricity])/size(temp2,1);
    
    
coor=round(cat(1,temporal.Centroid));
   

    %plot(coor(:,1),coor(:,2),'.');
   