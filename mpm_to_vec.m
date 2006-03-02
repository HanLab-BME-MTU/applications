clear
%mpm_variable = load('mpm_org.mat');
%mpm_variable = load('mpm_frommlg');
mpm_variable = load('mpm_yx_Feb24');


a=char(fieldnames(mpm_variable));
MPM = mpm_variable.(a);
clear mpm_variable 

img_h = 779;

% trace MPM
t_n =1;
track_length(1)=0;

for t=1:2:size(MPM,2)-2
    nr=1;
    for i=1:size(MPM,1)
        if MPM(i,t) ~= 0 & MPM(i,t+2) ~= 0
            vector(nr,:) = MPM(i,t:t+3);
            
%             vector(nr,1) = img_h - MPM(i,t+1);
%             vector(nr,2) = MPM(i,t);
%             vector(nr,3) = img_h - MPM(i,t+3);
%             vector(nr,4) = MPM(i,t+2);
            nr=nr+1;
        end
    end
    % save 
    num_ext = num2str(int8(t/2), '%03i')
    save(['vector' num_ext],'vector');
end
