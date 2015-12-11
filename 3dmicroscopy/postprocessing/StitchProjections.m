%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       small program to stitch axial MIPs 
%	manual selection of shift vector, averaging in overlap
%       primitive version with visual stitching for light-sheet imaging.
%       by Reto Fiolka  12-04-2015
%       
%   MIP are 2D images
%   shiftvec is a matrix that contains the lateral and vertical shifts
%   between MIP1 and Mip2 (first row in matrix) and Mip2 and Mip3 (second
%   row)
%                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MIP=StitchProjections(Mip1,Mip2,Mip3,shiftvec)



[zdim,xdim]=size(Mip1);


%%

Mip1=[zeros(2*zdim,xdim);Mip1];Mip2=[zeros(zdim,xdim);Mip2;zeros(zdim,xdim)];Mip3=[Mip3;zeros(2*zdim,xdim)];

%%
Mip1=Mip1-min(Mip1(:));Mip2=Mip2-min(Mip2(:));Mip3=Mip3-min(Mip3(:));

% manually set shift vectors

 xshift1=shiftvec(1,1); zshift1=shiftvec(1,2);
 xshift3=shiftvec(2,1); zshift3=shiftvec(2,2);

% xshift1=31;zshift1=-7; 
% xshift3=-7;zshift3=11;


Mip1s=circshift(Mip1,[zshift1,xshift1]);
scale1=mean2(Mip2(2*zdim+1+zshift1:2*zdim,:))/mean2(Mip1s(2*zdim+1+zshift1:2*zdim,:));
Mip1=Mip1*scale1;

Mip3s=circshift(Mip3,[zshift3,xshift3]);
scale3=mean2(Mip2(zdim+1:zdim+zshift3,:))/mean2(Mip3s(zdim+1:zdim+zshift3,:));
Mip3=Mip3*scale3;

MIP=circshift(Mip1,[zshift1,xshift1])+Mip2+circshift(Mip3,[zshift3,xshift3]);

MIP(zdim+1:zdim+zshift3,:)=MIP(zdim+1:zdim+zshift3,:)/2;
MIP(2*zdim+1+zshift1:2*zdim,:)=MIP(2*zdim+1+zshift1:2*zdim,:)/2;


MIP(2*zdim+1+zshift1,:)=Mip2(2*zdim+1+zshift1,:);









