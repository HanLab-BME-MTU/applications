clear
load('mpm_org.mat');

% trace MPM
t_n =1;
track_length(1)=0;
for i=1:size(MPM,1)
    j=1;
    while j < size(MPM,2)
        j_start =j;
        while j < size(MPM,2) && MPM(i,j) ~= 0 
            track_length(t_n) = track_length(t_n) + 1;
            j=j+1;
        end
        if track_length(t_n) < 200
            for j_del = j_start : j-1
                MPM(i,j_del) = 0;
                
                % clean M as well
                if mod(j_del,2)==0
                    M(i,1:4,(j_del)/2) = 0;
                end
                %elseif mod(j_del,2)~=0
                %    M(i,1:4,(j_del-1)/2) = 0;
                %end
            end
            if mod(j_del,2)==0 
                M(i,1:4,(j_del+2)/2) = 0;
            end

            t_n = t_n -1;
            
            
        end
        t_n = t_n + 1;
        track_length(t_n) = 0;
        j=j+2;
    end
    j=1;
end

t_n
[hi x_bin] = hist(track_length,30);
figure
plot(x_bin,hi);
max(track_length)
save('mpm_mod','M','MPM');
