function [coordnew,cc]= templfindertrack(coord,oldpic,newpic,levback,tempelcellmarker,newcell,newcelltempel,newcelltempelmarker,percentbackground,sizetemple,box_size_img)
% templfindertrack finds cells in a subsequent frame using template
%                  matching (correlation)
%
% SYNOPSIS       [coordnew,cc]= templfindertrack(coord,oldpic,newpic,levback,tempelcellmarker,newcell,newcelltempel,
%                                                newcelltempelmarker,percentbackground,sizetemple,box_size_img)
%
% INPUT 		 coord : coordinates of the cell we wish to find in the new image
% 				 oldpic : the image in which that cell is 
% 				 newpic : the image in which we wish to find the cell 
% 				 levback : approximated level of the background
% 				 tempelcellmarker : marker for old cells found by templates
% 				 newcell : identifier of new cells 
% 				 newcelltempel : identifier of new cells found by template
% 				 newcelltempelmarker : marker for new cells found by templates
% 			   	 percentbackground : threshold for ignoring areas within template and image (1+-percentbackground * levback = thresh)
% 				 sizetemple : the size we want the template to have
% 			 	 box_size_img : the size of the searcharea
% 
% OUTPUT         coordnew : found coordinates
%                cc : maximum of correlation
%
% DEPENDENCIES   templfindertrack uses {nothing}
%                                  
%                templfindertrack is used by { trackmater }
%                                   
% Colin Glass, Feb 04         

%make a template

[imgold_h,imgold_w]=size(oldpic);


xmin=round(coord(1,1)-(sizetemple-1)/2);
xmax=round(coord(1,1)+(sizetemple-1)/2);
ymin=round(coord(1,2)-(sizetemple-1)/2);
ymax=round(coord(1,2)+(sizetemple-1)/2);


if xmin < 1
      xmin =1;
end

if xmax > imgold_w;
      xmax = imgold_w;
end

if ymin <1
      ymin =1;
end

if ymax > imgold_h
      ymax= imgold_h;
end
   
%extract the pattern
template=oldpic(ymin:ymax,xmin:xmax);
 



%now we replace the background with random noise, so that we don't correlate
%background to background

barowtem=[];
bacolltem=[];

up= 1+percentbackground;
down= 1-percentbackground;

[barowtem,bacolltem]=find(template<(levback*up)&template>(levback*down));
   
for f=1:length(barowtem)
       template(barowtem(f),bacolltem(f))=rand(1)*255;
end
clear barowtem;
clear bacolltem;

[img_h,img_w]=size(newpic);




    
%size of the operation box, must be odd value!


%     if max(sizetemple)>box_size_img
%        box_size_img=max(sizetemple)+30;
%        if rem(box_size_img,2) == 0
%           box_size_img=box_size_img+1;
%        end
%     end
  
%get the blocks along the edge

%determine the length of the edge [pixel]
%[edge_l n]=size(edge);

x_img_1= round(coord(1,1)-(box_size_img-1)/2);
x_img_2= round(coord(1,1)+(box_size_img-1)/2);
y_img_1= round(coord(1,2)-(box_size_img-1)/2);
y_img_2= round(coord(1,2)+(box_size_img-1)/2);
   
if x_img_1 < 1
      x_img_2 = x_img_2 - x_img_1 + 1;
      x_img_1 = 1;
end

if x_img_2 > img_w
      x_img_1 = x_img_1 - x_img_2 + img_w;
      x_img_2 = img_w;    
end

if y_img_1 <1
      y_img_2 = y_img_2 - y_img_1 + 1;
      y_img_1 = 1;
end

if y_img_2 > img_h
      y_img_1 = y_img_1 - y_img_2 + img_h;
      y_img_2 = img_h;
end
   
   
   
%extract part of the image
searcharea=newpic(y_img_1:y_img_2, x_img_1:x_img_2);

barow=[];
bacoll=[];
[barow,bacoll]=find(searcharea<(levback*1.2)&searcharea>(levback*0.8));
   
for f=1:length(barow)
       searcharea(barow(f),bacoll(f))=rand(1);
end
clear barow;
clear bacoll;

%do a normalized crosscorrelation
img_corr=[];
sizecorr=[];

img_corr=normxcorr2(template,searcharea);

sizecorr=size(img_corr);

%now we superimpose a distance criteria. First we have to build a matrix of
%the same size as img_corr, then adjust it's values according to the
%distance from where the cell was in the last picture
dist = zeros(sizecorr(1,1),sizecorr(1,2));

for f = 1:sizecorr(1,1)
       for g = 1:sizecorr(1,2)
               
               %the following means: if the distance to the initial
               %coordinates is smaller than 20, the value is one
               if (sqrt((g+x_img_1-round((sizetemple-1)/2)-coord(1,1))^2+...
                       (f+y_img_1-round((sizetemple-1)/2)-coord(1,2))^2))<20
                        dist(f,g)=1;
                        
               %if it is bigger than 20, the value slowly decreases         
               else dist(f,g) = 1-((sqrt((g+x_img_1-round((sizetemple-1)/2)-coord(1,1))^2+...
                                (f+y_img_1-round((sizetemple-1)/2)-coord(1,2))^2)-20)/80);
               end
              
               %negative values become zero
               if dist(f,g)<0
                       dist(f,g)=0;
               end
       end
end

%multiply img_corr with dist, thus having intrduced a distance criteria
img_cor = img_corr.*dist;

clear dist;

% % % %mask the central peak
% % % %[corr_h,corr_w,l]=size(img_corr);
% % % %img_corr((corr_h-1)/2-3:(corr_h-1)/2+3,(corr_w-1)/2-3:(corr_w-1)/2+3)=0;
  

%find the coordinates of the absolute intensity maxima in the image
[c xi]=max(img_cor);
%%%%%werte(i)=max(c);
[cc xii]=max(c);
yjj=xi(xii);
   

   
   
%add together everything that's necessary to receive the coordinates of the
%found spot in the frame (not within img_corr)
coordnew(1,1)=xii + x_img_1 - round(size(template,1)/2);      %%%+round(difftemplecentcentroid(i,1));
coordnew(1,2)=yjj + y_img_1 - round(size(template,2)/2);   %%%+round(difftemplecentcentroid(i,2));

if coordnew(1,1) < 1
      coordnew(1,1) = 1;
end

if coordnew(1,1) > (img_w-1)
      coordnew(1,1)= img_w-1;    
end

if coordnew(1,2) < 1
      coordnew(1,2) = 1;
end

if coordnew(1,2) > (img_h-1)
      coordnew(1,2) = img_h-1;
end
   

coordnew(1,1) = coordnew(1,1)+tempelcellmarker;
coordnew(1,2) = coordnew(1,2)+tempelcellmarker;


%mark cells that were new cells (or new cells propagated by 
%templates) as new cells propagated by templates
coordflo = floor(coord+0.1);
what = round(10*(coord-coordflo));

if what==newcell | what==newcelltempel
        coordnew = floor(coordnew)+newcelltempelmarker;
end
    
   
% % %    
% % %    %figure, imshow(img_corr), title('corralation')
% % %    %hold on;
% % %    %plot(xii,yjj,'.');
% % %    %hold off;
% % %    
% % %    %find the value of the second biggest maximum
% % %    %img_corr((yjj-20):(yjj+20),(xii-20):(xii+20))=0;
% % %   %[c2 xi2]=max(img_corr);
% % %   % werte2(i)=max(c2);
% % %  %  werte2(i)=0;
% % %    
% % %    %%[cc2 xii2]=max(c2);
% % %    %%yjj2=xi2(xii2);
% % %    %%figure, imshow(img_corr), title('corralation')
% % %    %%hold on;
% % %    %%plot(xii2,yjj2,'.');
% % %    %%hold off;
% % %    
% % %    %[cc xii]=max(c2);
% % %    %yjj=xi(xii);
% % %    %xuhm=xii-round((sizetemple(1,1)-1)/2);
% % %    %yuhm=yjj-round((sizetemple(1,2)-1)/2)
% % %    %figure, imshow(searcharea)
% % %   %hold on;
% % %   % plot(xuhm,yuhm,'.');
% % %    %hold off;
% % %   
% % % % % else
% % % % %     coordnew(1,1)=0;
% % % % %     coordnew(1,2)=0;
% % % % %     cc=0;
% % % % % end