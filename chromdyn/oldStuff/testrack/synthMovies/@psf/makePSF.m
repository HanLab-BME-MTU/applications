function psfdata=makePSF(ps,dim,samp,pos)
%Create a psf intensity matrix
xd=floor(dim(1)/2);
yd=floor(dim(2)/2);
zd=floor(dim(3)/2);
i=0;
psfdata=zeros(dim);
h = waitbar(0,'Please wait...');
for x=-xd:xd
   for y=-yd:yd
      for z=-zd:zd
         psfdata(x+xd+1,y+yd+1,z+zd+1)=psfValue(ps,samp(1)*(x-pos(1)),samp(1)*(y-pos(2)),samp(2)*(z-pos(3)));
      end;
      i=i+1;
      waitbar(i/(dim(1)*dim(2)),h)
   end;
end;
close(h);
   
