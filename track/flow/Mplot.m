function Mplot=Mplot(M)
 M=double(M);
 uzt1=size(M,1);
 nc112=size(find(M(:,3)==0));
 nc212=size(find(M(:,1)==0));
 uzt1=uzt1(1)-nc112-nc212;
 Mbase=M(1:uzt1,:);
 Ma=M(uzt1+1:uzt1+nc112,:); 
 Maa=M(uzt1+nc112+1:uzt1+nc112+nc212,:);
 
hold on
axis ij
QUIVER(Mbase(:,2),Mbase(:,1),Mbase(:,4)-Mbase(:,2),Mbase(:,3)-Mbase(:,1),0);

plot(Ma(:,2),Ma(:,1),'b*'); %plot the notarget
plot(Maa(:,4),Maa(:,3),'r*'); %plot the no source

