function plot3Data(img,lthres,uthres)
% PLOT3DATA plots points in 3d image stack
%
% SYNOPSIS plot3Data(img,lthres,uthres)
%
% INPUT img : image=[x,y,z]
%             lthres= lower threshold
%             upthres= (optional) upper threshold 

% 12/10/00 dT

%check for uthres
if nargin==2
   uthres=inf;
end;
wrnH=-1;

imC=find((img>lthres) & (img<uthres));
%create cube
vx=[0 1 1 0 0 1 1 0]';
vy=[0 0 1 1 0 0 1 1]';
vz=[0 0 0 0 1 1 1 1]';
vert=[vx vy vz];
face = [1 2 6 5; 2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4; 5 6 7 8];
[x,y,z]=ind2sub(size(img),imC);
img=img/max(img(:));
colordef black;

%limit the number of cubes to 5000
if length(x)>5000
   dumF=gcf;
   wrnH=warndlg({'Attempt to draw more than 5000 cubes.' 'Only 5000 are drawn!'},'plot3data');
   pause(0.1);
   figure(dumF);
   %resample
   ind=round(1:length(x)/5000:length(x));
   x=x(ind);
   y=y(ind);
   z=z(ind);
end;

%draw cubes
for l = 1:length(x)
   c=img(x(l),y(l),z(l));
   vert(:,1)=vx*c+x(l);
   vert(:,2)=vy*c+y(l);
   vert(:,3)=vz*c+z(l);
   h=patch('Vertices',vert,'Faces',face,'FaceVertexCData',ones(8,3)*c,...
      'FaceColor','interp','LineStyle','none');
end;
rotate3d on;
if ishandle(wrnH)
   close(wrnH);
end;
