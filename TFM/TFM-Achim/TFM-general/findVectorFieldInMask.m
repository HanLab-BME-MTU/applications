function [inpos,invec]=findVectorFieldInMask(pos,vec,bwMask)
% pos and bwMask must live in the same coordinate system. In particular,
% bwMask(1,1)=[1,1] the origin of the coordinate system.
size1= max(fliplr(pos));
size2=size(bwMask);
ImgSize=max(vertcat(size1,size2));
linIndicesFpos=sub2ind(ImgSize,pos(:,2),pos(:,1));
fx_mat=zeros(ImgSize);
fy_mat=zeros(ImgSize);

fx_mat(linIndicesFpos)=vec(:,1);
fy_mat(linIndicesFpos)=vec(:,2);

% extend the mask if smaller than the image:
check_mat1=bwMask;
if numel(check_mat1)<prod(ImgSize)
    check_mat1(ImgSize(1),ImgSize(2))=0;
end
%check_mat1= logical(cellMask);

% find those entries that are in both lists:
check_mat2= false(ImgSize);
check_mat2(linIndicesFpos)=true;
check_mat=(check_mat1 & check_mat2);

[cell_ypos,cell_xpos] = ind2sub(ImgSize,find(check_mat));
cell_xvec=fx_mat(check_mat);
cell_yvec=fy_mat(check_mat);

% the vectors located inside are:
inpos=horzcat(cell_xpos,cell_ypos);
invec=horzcat(cell_xvec,cell_yvec);