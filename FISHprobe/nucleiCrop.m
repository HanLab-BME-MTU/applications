function nucleus = nucleiCrop(num,channel,imInput,imLabelCellSeg)

nucMask=zeros(size(imInput));

% Given the segmented nuclei number
nucMask(imLabelCellSeg==num)=1;

% Modify the dilation factors as is
% se=strel('ball',10,3);
[a,b,c]=ndgrid(-10:10,-10:10,-3:3);
se=strel(sqrt(a.^2+b.^2+c.^2));
nucMask2=imdilate(nucMask,se);

nucData=channel;
nucData(nucMask2==0)=0;

i=1;
j=1;
k=1;
while i<=size(nucData,1)
    if min(min(nucData(i,:,:)==0))==1
        nucData(i,:,:)=[];
        i=i-1;
    end
    i=i+1;
end

while j<=size(nucData,2)
    if min(min(nucData(:,j,:)==0))==1
        nucData(:,j,:)=[];
        j=j-1;
    end
    j=j+1;
end

while k<=size(nucData,3)
    if min(min(nucData(:,:,k)==0))==1
        nucData(:,:,k)=[];
        k=k-1;
    end
    k=k+1;
end

nucleus=nucData;
% 3D Plot

v=nucData;
if size(v,1) < size(v,2)
    v(size(v,1)+1:size(v,2),:,:)=0;
else if size(v,1) > size(v,2)
        v(:,size(v,2)+1:size(v,1),:)=0;
    end
end

[x,y,z]=meshgrid(1:size(v,1),1:size(v,2),1:size(v,3));

figure,
p = patch(isosurface(x,y,z,v));
% why not -3? see doc isosurface
isonormals(x,y,z,v,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,39/152])
view(3); axis tight
camlight
lighting gouraud
title(strcat('Nucleus',num2str(num)));
% set white background
set(gcf,'color','white')