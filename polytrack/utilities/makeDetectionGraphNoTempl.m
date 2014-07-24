cd ('/home/kerstens/journal_paper/tracking_quality/5min/results_start1_s12_t1/info');

% Init matrices
pass1 = zeros(193,1);
pass2 = zeros(193,1);
pass3 = zeros(193,1);
pass4 = zeros(193,1);
pass5 = zeros(193,1);

for i = 1 : 193
    
    % Cells found by kmeans method
    error = 0;
    inputFile = ['nucleiCoord_pass1_' num2str(i)];
    try
        load (inputFile);
    catch
        disp(['Cannot load mat file ' inputFile '.']);
        error = 1;
    end
    if ~error
        pass1(i) = length(nucCoord);
    end
    
    % Cells found by looking at halos
    error = 0;
    inputFile = ['nucleiCoord_pass2_' num2str(i)];
    try
        load (inputFile);
    catch
        disp(['Cannot load mat file ' inputFile '.']);
        error = 1;
    end
    if ~error
        pass2(i) = length(newCoord);
    end

    % Cells found by cluster size information
    error = 0;
    inputFile = ['nucleiCoord_pass3_' num2str(i)];
    try
        load (inputFile);
    catch
        disp(['Cannot load mat file ' inputFile '.']);
        error = 1;
    end
    if ~error
        pass3(i) = length(newCoord);
    end
end

% Load the MPM and get final coords
load ('MPMBeforeProcessing.mat');
for i = 1 : 193
    pass4(i) = length(find(MPM(:,i*2)));
end

% Load the MPM and get coords after min track linking (2 points min.)
load ('../MPM.mat');
for i = 1 : 193
    pass5(i) = length(find(MPM(:,i*2)));
end

% Read in the manually counted cells
if exist('121404_s12_manual_count.csv')
    manCount = csvread('121404_s12_manual_count.csv');
    manCount = manCount(:,2);
end

% Open a figure for the cell amounts
f1 = figure; plot(pass1,'r');
xlabel ('Frames');
ylabel ('Nr of cells');
hold on;
plot(pass2,'g');
plot(pass3,'b');
plot(pass4,'k');
plot(pass5,'c');
if exist('manCount','var')
    plot(manCount,'m');
    legend('Kmeans clust.','Halos','Cluster size','Gap closing',...
           'After linking','Manual count','Location','NorthWest');
else
    legend('Kmeans clust.','Halos','Cluster size','Gap closing',...
           'After linking','Location','NorthWest');
%          'Min. track length','Location','NorthWest'); 
end
hold off;

% Put a title
title('Detection of cells');


% Open a figure for the percentage of cells found
if exist('manCount','var')
    percentage = (pass4./manCount).*100;
    indx = find(percentage>100);
    percentage(indx) = 100;
    f1 = figure; plot(percentage,'r');
    xlabel ('Frames');
    ylabel ('%');
    axis([1 193 0 110]);
end

% Put a title
title('Percentage of cells found');


% % Open a figure for the delta amounts
% pass12 = pass2-pass1;
% pass23 = pass3-pass2;
% pass34 = pass4-pass3;
% %pass45 = pass5-pass4;
% f2 = figure; plot(pass12,'r');
% hold on;
% plot(pass23,'g');
% plot(pass34,'b');
% %plot(pass45,'c');
% legend('Kmeans clust.','Halos','Cluster size','Gap closing','Location','NorthWest');
% hold off;
% 
% % Put a title
% title('Delta amounts');


% Open a figure for the division amounts
pass14 = pass1./pass4;
pass24 = (pass2-pass1)./pass4;
pass34 = (pass3-pass2)./pass4;
pass4e = (pass4-pass3)./pass4;
f3 = figure; plot(pass14.*100,'r');
xlabel ('Frames');
ylabel ('%');
hold on;
plot(pass24.*100,'g');
plot(pass34.*100,'b');
plot(pass4e.*100,'k');
legend('Perc. kmeans clust.','Perc. from halos','Perc. from cluster size',...
       'Perc. from gap closing');
hold off;

% Put a title
title('Divided amounts');
