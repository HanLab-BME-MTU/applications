function Y2=extendD(Y,dim)
szmin=min([size(Y)'; dim']);
Y2=zeros(dim);
Y2(1:szmin(1),1:szmin(2),1:szmin(3))=Y(1:szmin(1),1:szmin(2),1:szmin(3));

% Mirror boundary conditions
if size(Y,1)<dim(1),
  sizediff=dim(1)-size(Y,1)-1;
  Y2(size(Y,1)+1:dim(1),:,:)=Y2(size(Y,1):-1:size(Y,1)-sizediff,:,:);
end;
if size(Y,2)<dim(2),
  sizediff=dim(2)-size(Y,2)-1;
  Y2(:,size(Y,2)+1:dim(2),:)=Y2(:,size(Y,2):-1:size(Y,2)-sizediff,:);
end;
if size(Y,3)<dim(3),
  sizediff=dim(3)-size(Y,3)-1;
  Y2(:,:,size(Y,3)+1:dim(3))=Y2(:,:,size(Y,3):-1:size(Y,3)-sizediff);
end;