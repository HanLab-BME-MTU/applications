%%
clear all
tic
xdim=512;zdim=126;
begin=1

PSF=zeros(zdim,xdim,xdim);
for i=1:zdim; 
   temp=double(imread('/home2/rfiolka/200nmbeads/Dense100micronRange200nmStep/PSF.tif',i));
   PSF(i,:,:)=temp(end-511:end,end-511:end);

end

%%
OTF=(fftshift(fftn(PSF)));



%% rotational averaging
OTFtemp=zeros(zdim,xdim-1,xdim-1);
OTFrot=zeros(zdim,xdim,xdim);

for k=1:zdim
test=squeeze(OTF(k,:,:));


len=xdim/2-1;c1=xdim/2+1;,c2=xdim/2+1;
  
a2sub = abs(test((c2-len):(c2+len),(c1-len):(c1+len)));

dist = ((-len:len)'*(ones(1,2*len+1))).^2;dist = (dist+dist').^.5;
dist = round(dist)+1; 

ravg = 1:len+1;
for i=1:len+1
    ravg=mean((a2sub(dist==i)));
a2sub(dist==i)=ravg;
end

OTFtemp(k,:,:)=a2sub;
end

OTFrot(:,2:end,2:end)=OTFtemp;

toc
%%
save(['OTFrot' num2str(zdim) 'slicesZ' num2str(xdim)  'compress.mat'],'OTFrot','-v7.3')



%%






