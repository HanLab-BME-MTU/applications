% test of plotting
% Liya Ding

close all;

% data generate

t=0:0.2:3*pi;

Mean_x = 5*sin(t);

error_x = randn(10, numel(Mean_x));

Data_Input = error_x + repmat(Mean_x, [10,1]);

%% Calculate mean, std, median, percentiles

Data_Mean = mean(Data_Input,1);

Data_std = (std(Data_Input));

Data_Median = median(Data_Input);

Data_1q = (quantile(Data_Input,0.25))';
Data_3q = (quantile(Data_Input,0.75))';

Data_10p = (quantile(Data_Input,0.10))';
Data_90p = (quantile(Data_Input,0.90))';


%% Ploting

% Data
figure(1); plot(Data_Input');


% Mean with 2 times STD
figure(2);
plot(Data_Mean,'LineWidth',6,'color','b');

hold on;

plot(Data_Mean-2*Data_std,'LineWidth',2,'color','m');
plot(Data_Mean+2*Data_std,'LineWidth',2,'color','m');


% Mean with 2 STD with blocks
figure(3);

for i = 1 :  numel(Data_Mean)
    
    X = [i-0.5 i-0.5 i+0.5 i+0.5 i-0.5];
    Y = [Data_Mean(i)-2*Data_std(i) Data_Mean(i)+2*Data_std(i) ...
        Data_Mean(i)+2*Data_std(i) Data_Mean(i)-2*Data_std(i) ...
        Data_Mean(i)-2*Data_std(i)];   
    
    patch(X,Y,'m');
    hold on;
end

plot(Data_Mean,'LineWidth',6,'color','b');


% Median with percentiles 10%, 90%
figure(4);

hold on;

plot(Data_10p,'LineWidth',2,'color','m');
plot(Data_90p,'LineWidth',2,'color','m');

plot(Data_Median,'LineWidth',6,'color','b');


% Median with quaters, blocks
figure(5);

hold on;

for i = 1 :  numel(Data_Mean)
    
    X = [i-0.5 i-0.5 i+0.5 i+0.5 i-0.5];
    Y = [Data_1q(i) Data_3q(i) ...
        Data_3q(i) Data_1q(i) ...
        Data_1q(i)];
        
    patch(X,Y,'m');
    
end

plot(Data_Median,'LineWidth',6,'color','b');

% Median with percentiles, blocks
figure(6);

hold on;

for i = 1 :  numel(Data_Mean)
    
    X = [i-0.5 i-0.5 i+0.5 i+0.5 i-0.5];
    Y = [Data_10p(i) Data_90p(i) ...
        Data_90p(i) Data_10p(i) ...
        Data_10p(i)];
        
    patch(X,Y, [0.7 0.7 0.99]);
    
end

plot(Data_Median,'LineWidth',6,'color','b');

