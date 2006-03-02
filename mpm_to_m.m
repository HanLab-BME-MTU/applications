clear
%mpm_variable = load('mpm_org.mat');
mpm_variable = load('mpm_yx_Feb24');

a=char(fieldnames(mpm_variable));
MPM = mpm_variable.(a);
clear mpm_variable 

for i=1:size(MPM,1)
    nr=1;   
    for t=1:2:size(MPM,2)-2
    	M(i,:,nr) = MPM(i,t:t+3);
        nr=nr+1;
    end
end
save('mpm_yx_Feb24_mod','M','MPM');