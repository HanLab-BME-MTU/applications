function [coor,halo]= halosfind(info,erodedisksize,halolevel)




index=[];
halo=zeros(size(info,1),size(info,2));
img_labels=[];
haloshlabel=[];
theones=[];

%look for pixels that are above a certain value
index=find(info >halolevel);
halo(index)=1;
clear index;

%close small gaps
se = strel('disk',3);
halo=imclose(halo,se);

%label groups
img_labels= bwlabel(halo);


%erode big time. In this way even big groups are reduced to nothing, if
%their area is not coherent
Se = strel('disk',erodedisksize);
haloshlabel= imerode(img_labels,Se);

clear img_labels;

%haloshrunk= imdilate(haloshrunk,Se);
%info is the picture, lev was determined from values taken from the picture

% % se = strel('disk',4);       % disk, radius 4
% % 
% % haloshrunk= imdilate(haloshr,se);

%figure, imshow(haloshlabel);

%see what groups are still present
theones=unique(haloshlabel);
theones(1,:)=[];



% % % % % % % % % huntdown=regionprops(img_labels,'Area','PixelList');%%%%%,'MajorAxisLength','Centroid','Eccentricity');
% % % % % % % % % %look up regionprops to add qualities of interest
% % % % % % % % % 
% % % % % % % % % %kick out really small and big areas
% % % % % % % % % img_area=[huntdown.Area];
% % % % % % % % % destroy=find(or(img_area < 270,img_area > 6000));
% % % % % % % % % temp=huntdown(destroy);
% % % % % % % % % 
% % % % % % % % % destcoord=cat(1,temp(:).PixelList);
% % % % % % % % % for o=1:size(destcoord,1)
% % % % % % % % % regmax(destcoord(o,2),destcoord(o,1))=0;
% % % % % % % % % end
% % % % % % 
% % % % % % 
% % % % % % %se = strel('disk',4);       % disk, radius 4
% % % % % % 
% % % % % % 
% % % % % % %figure, imshow(dilregmax);
% % % % % % 
% % % % % % 
% % % % % % %img_labels= bwlabel(dilregmax);


%calculate properties
realnucloi=regionprops(haloshlabel,'Area','PixelList','MajorAxisLength','Centroid','Eccentricity');
    %look up regionprops to add properteties of interest

%get the centroids. Note that the centroids are of the eroded picture.
%Generally properties of the eroded picture are NOT representative!
yeah=realnucloi(theones);
coor=round(cat(1,yeah.Centroid));


% % % % % % % % % % % % % % %kick out really small and big areas
% % % % % % % % % % % % % % img_areare=[realnucloi.Area];
% % % % % % % % % % % % % % eccent=[realnucloi.Eccentricity];
% % % % % % % % % % % % % % relevant=find(img_areare > 300 & eccent<0.7);
% % % % % % % % % % % % % % %relevant2=find(eccent>0.6
% % % % % % % % % % % % % % %relevant=
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % tempor=realnucloi(relevant);
% % % % % % % % % % % % % % clear img_areare;
% % % % % % % % % % % % % % img_areare=[tempor.Area];
% % % % % % % % % % % % % % relevant2=find(img_areare < 3000);
% % % % % % % % % % % % % % temporal=tempor(relevant2);
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % clear img_area;
% % % % % % % % % % % % % % % % img_area=[temp.Area];
% % % % % % % % % % % % % % % % relevant2=find(img_area < 3000);
% % % % % % % % % % % % % % % % temp2=temp(relevant2);
% % % % % % % % % % % % % % % % %meanecc=sum([temp2.Eccentricity])/size(temp2,1);
% % % % % % % % % % % % % %     
% % % % % % % % % % % % % %     
% % % % % % % % % % % % % % coor=round(cat(1,temporal.Centroid));
% % % % % % % % % % % % % %    
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % %     %plot(coor(:,1),coor(:,2),'.');
% % % % % % % % % % % % % %    