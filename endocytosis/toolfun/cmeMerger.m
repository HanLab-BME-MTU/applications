% This function provides the factor for merging data on different days.
% The factor is an array of numbers for normalizing the maximum intensity
% thresholds, which are required for cmeAnalysis. 
% Once the factor, and the MaxIntThreshold of day 1 is obtained, then the
% threshold on different days can be calculated using the factor:
% MaxIntThreshold(2) = MaxIntThreshold(1) * factor(1) / factor(2), for day2
% MaxIntThreshold(3) = MaxIntThreshold(1) * factor(1) / factor(3), for day3
% and so forth
% 1. Call the function: yyy = cmeMerger(xxx), xxx is the number of days to be merged, yyy is the factor. 
% 2. Choose preferred output folder
% 3. factor will be stored as 'factor.mat'
%=============================================================
% Optionally, 'numMinErrStep' is the step number to minimize the difference of the cumulative
% distribution functions among days
%=============================================================
% Author: Xinxin Wang, Danuser Lab

function factor = cmeMerger(num_condition, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('num_condition', @(x) isnumeric(x));
ip.addOptional('MaximumIntensity', 500, @(x) isnumeric(x));
ip.addOptional('numMinErrStep', 100, @(x) isnumeric(x));
ip.addOptional('resultDir', [], @(x) ischar(x));
ip.parse(num_condition, varargin{:});


data_master_dir = ip.Results.resultDir;
n_condition = ip.Results.num_condition;
h = 1.0;
n_h = ip.Results.MaximumIntensity / h;
n_opt = ip.Results.numMinErrStep;

A_min = 0;

cpdf = zeros(n_condition, n_h);
cpdf_ave = zeros(1, n_h);
diff_tot_temp = zeros(1, n_h);
diff_tot = 1000.*ones(1, n_h);
diff = zeros(n_condition, n_h);
factor = ones(n_condition, 1);
x_ave = (h:n_h);
x = zeros(n_condition, n_h);
for i = 1:n_h
x(i,:) = (h:n_h);
end
%============================================================
if isempty(data_master_dir)
    data_master_dir = uigetdir(pwd, 'Select output folder:');
    if data_master_dir==0
        return
    else
        data_master_dir = [data_master_dir filesep];
    end
end
if ~strcmp(data_master_dir(end), filesep)
    data_master_dir = [data_master_dir filesep];
end
cd(data_master_dir);
%============================================================
num_movie = zeros(1, n_condition);
for i_condition = 1:n_condition
    disp('Load the data for day #');
    disp(i_condition);
    data = loadConditionData();  
    data_size = size(data);
    num_movie(i_condition) = data_size(2);
    if i_condition == 1
        data_tot = data;
    else
        data_tot = poolExperiments(data_tot, data);
    end
end
cd(data_master_dir);
%save('merging_result.mat','data_tot','num_movie');
%============================================================
for i_condition = 1:n_condition
%============================================================
if i_condition == 1
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

for i_data = 1:n_data(2)
 disp('movie#:');
 disp(i_data);
 track = loadTracks(data(i_data));
 n_track = size(track);
 for i = 1:n_track(2)
     for j = 1: n_h
         bound = A_min+j*h;
         n_frame = size(track(i).A);
         for k = 1:  n_frame(2)
             if ((track(i).A(k) < bound) && (track(i).A(k) > A_min))
                 cpdf(i_condition,j) = cpdf(i_condition,j) +1 ;
             end
         end
     end
 end
end
end
%==============================================================
for i_condition = 1:n_condition
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
for i_opt = 1:n_opt  
    disp('Percentage finished for calculation:');
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
