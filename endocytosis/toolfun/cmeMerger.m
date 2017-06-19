

h = 1.0;
n_h = 500;
A_min = 0;
A_max = A_min + n_h*h;
prompt = 'How many days of data do you need to merge? ';
n_condition = input(prompt);
cpdf = zeros(n_condition, n_h);
cpdf_ave = zeros(1, n_h);
diff_tot_temp = zeros(1, n_h);
diff_tot = 1000.*ones(1, n_h);
diff = zeros(n_condition, n_h);
factor = ones(n_condition, 1);
factor_temp = ones(n_condition, 1);
x_ave = (h:n_h);
x = zeros(n_condition, n_h);
for i = 1:n_h;
x(i,:) = (h:n_h);
end
%============================================================
prompt = 'Enter your prefered folder name:';
master_folder_name = input(prompt,'s');
mkdir(master_folder_name);
cd(master_folder_name);
data_master_dir = pwd;
%============================================================
num_movie = zeros(1, n_condition);
for i_condition = 1:n_condition;
    disp('Load the data for day #');
    disp(i_condition);
    data = loadConditionData();  
    data_size = size(data);
    num_movie(i_condition) = data_size(2);
    if i_condition == 1;
        data_tot = data;
    else
        data_tot = poolExperiments(data_tot, data);
    end
end
cd(data_master_dir);
save('data_tot.mat','data_tot','num_movie');
%============================================================  
%data_dir = cell(1,n_condition);
%for i_condition = 1:n_condition;
%    if i_condition == 1;
%     data_tot_start = 1;
%    else
%     data_tot_start = sum(num_movie(1:(i_condition-1)))+1;
%    end
%data = data_tot(data_tot_start);
%data_dir{i_condition} = data_tot(i_condition).source;
%end
%============================================================
cd(data_master_dir);
for i_condition = 1:n_condition;
    mkdir(num2str(i_condition));
end
output_dir = cell(1,n_condition);
for i_condition = 1:n_condition;
    output_dir{i_condition} = [data_master_dir filesep num2str(i_condition)];
end
%============================================================
for i_condition = 1:n_condition;
    cd(output_dir{i_condition});
    for i_movie = 1:num_movie(i_condition);
        mkdir(num2str(i_movie));
    end
end
%============================================================
%============================================================
for i_condition = 1:n_condition;
%============================================================
if i_condition == 1;
data_tot_start = 1;
else
data_tot_start = sum(num_movie(1:(i_condition-1)))+1;
end
data_tot_end = sum(num_movie(1:(i_condition)));
data = data_tot(data_tot_start:data_tot_end);
disp('reading condition#:');
disp(i_condition);
n_data = size(data);
%============================================================

for i_data = 1:n_data(2);
 disp('movie#:');
 disp(i_data);
 track = loadTracks(data(i_data));
 n_track = size(track);
 for i = 1:n_track(2);
     for j = 1: n_h;
         bound = A_min+j*h;
         n_frame = size(track(i).A);
         for k = 1:  n_frame(2);
             if ((track(i).A(k) < bound) && (track(i).A(k) > A_min))
                 cpdf(i_condition,j) = cpdf(i_condition,j) +1 ;
             end
         end
     end
 end
end
end
%==============================================================
for i_condition = 1:n_condition;
cpdf(i_condition,:) = cpdf(i_condition,:)/cpdf(i_condition,n_h);
end
%==============================================================
for j = 1:n_h
    for i = 1:n_condition
        cpdf_ave(j) = cpdf_ave(j) + cpdf(i,j);
    end
    cpdf_ave(j) = cpdf_ave(j)/n_condition;
end
plot ( x(1,:), cpdf(1,:), 'r-',  x(2,:), cpdf(2,:), 'g-', x_ave, cpdf_ave(:), 'b-', 'Linewidth', 2 );
%==============================================================
n_opt = 100;
for i_opt = 1:100  
    disp('Percentage finished for analysis:');
    disp(i_opt/n_opt*100);
    %-------------------------------------------------------
            if i_opt == 1
				factor_temp = ones(n_condition, 1);
            else
                factor_temp = ones(n_condition, 1) + 0.1*normrnd(0,1);
            end
    %-------------------------------------------------------
    for i = 1:n_condition
        for j = 2:n_h-1
            diff_grid_min = 10.;
            for k = 1:n_h
                diff_grid = abs(x_ave(j) - x(i,k) * factor_temp(i));				
			    if diff_grid < diff_grid_min
                   diff_grid_min = diff_grid;
                   if x_ave(j) - x(i,k) * factor_temp(i) < 0.
                       cpdf_ave_new = cpdf_ave(j) + (cpdf_ave(j + 1) - cpdf_ave(j)) / (x_ave(j + 1) - x_ave(j)) * (x(i,k) * factor_temp(i) - x_ave(j));
                   else
                       cpdf_ave_new = cpdf_ave(j) + (cpdf_ave(j - 1) - cpdf_ave(j)) / (x_ave(j - 1) - x_ave(j)) * (x(i,k) * factor_temp(i) - x_ave(j));
                   end
                   diff(i,j) = (cpdf_ave_new - cpdf(i,k));
                end
            end
        end
    end
    %---------------------------------------------------------------
    for i = 1:n_condition
        diff_tot_temp(i) = 0.;
        for j = 1:n_h       
            diff_tot_temp(i) = diff_tot_temp(i) + abs(diff(i,j));
        end
        		if diff_tot_temp(i) < diff_tot(i)
			     diff_tot(i) = diff_tot_temp(i);
			     factor(i) = factor_temp(i);
                end
    end
   %---------------------------------------------------------------
end
%==============================================================================
cd(data_master_dir);
save('factor.mat','factor');
%==============================================================================
figure();
    	for i = 1:n_condition
		for k = 1:n_h
			x(i,k) = x(i,k)* factor(i);
        end
        end
plot ( x(1,:), cpdf(1,:), 'r-',  x(2,:), cpdf(2,:), 'g-', x_ave, cpdf_ave(:), 'b-', 'Linewidth', 2 );
%==============================================================================
