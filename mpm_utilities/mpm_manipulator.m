clear
mpm_variable = load('mpm.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine if ther is MPM and M variable
% stored or only MPM
a=char(fieldnames(mpm_variable));
if size(a,1) > 1
    MPM = mpm_variable.(a(2,:));
else
    MPM = mpm_variable.(a);
end
clear mpm_variable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trace MPM
t_n =1;

t_n_d = 1;
track_length(1)=0;
track_data=zeros(200000,3);

for i=1:size(MPM,1)
    i
    j=1;
    while j < size(MPM,2)
        j_start =j;
        
        % save start point of track
        track_data(t_n,1) = MPM(i,j);
        track_data(t_n,2) = MPM(i,j+1);
        
        while j < size(MPM,2) && MPM(i,j) ~= 0 
            track_length(t_n) = track_length(t_n) + 1;
            j=j+1;
        end
        
        % save length of track
        track_data(t_n_d,3) = track_length(t_n);
        t_n_d = t_n_d+1;
        
        if 0
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
        end % if 0
        t_n = t_n + 1;
        track_length(t_n) = 0;
        j=j+2;
    end
    j=1;
end

% plot track lenght map
track_map = zeros(max(track_data(:,1)),max(track_data(:,2)));
track_map_s = zeros(max(track_data(:,1)),max(track_data(:,2)));
max_length = max(track_data(:,3));
%cmap = colormap(jet(max_length));
cmap = colormap(jet(100));
for i=1:size(track_data,1)
     track_map(track_data(i,1),track_data(i,2),1) = cmap(track_data(i,3),1);
     track_map(track_data(i,1),track_data(i,2),2) = cmap(track_data(i,3),2);
     track_map(track_data(i,1),track_data(i,2),3) = cmap(track_data(i,3),3);
     %track_map(track_data(i,1),track_data(i,2),1) = 1;
     %track_map(track_data(i,1),track_data(i,2),2) = 1;
     %track_map(track_data(i,1),track_data(i,2),3) = 1;
end
figure,imshow(track_map,[]);


track_data_short = track_data(track_data(:,3) <= 3,:);
track_data_long  = track_data(track_data(:,3) > 3,:);
for i=1:size(track_data_short,1)
     track_map_s(track_data_short(i,1),track_data_short(i,2),1) = 1;
     track_map_s(track_data_short(i,1),track_data_short(i,2),2) = 0;
     track_map_s(track_data_short(i,1),track_data_short(i,2),3) = 0;
end
for i=1:size(track_data_long,1)
     track_map_s(track_data_long(i,1),track_data_long(i,2),1) = 1;
     track_map_s(track_data_long(i,1),track_data_long(i,2),2) = 0;
     track_map_s(track_data_long(i,1),track_data_long(i,2),3) = 1;
end
figure,imshow(track_map_s,[]);


if 0
    t_n
    [hi x_bin] = hist(track_length,30);
    figure
    plot(x_bin,hi);
    max(track_length)
    save('mpm_mod','M','MPM');
end




