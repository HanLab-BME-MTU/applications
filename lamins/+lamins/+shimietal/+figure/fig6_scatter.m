%% LB1null vs WT, Edges per face vs Area
h = figure;
scatter(stats(20).area(1).all.data,stats(20).edgesPerFace(1).all.data,'Marker','o')
hold on;
scatter(stats(7).area(1).all.data,stats(7).edgesPerFace(1).all.data,'Marker','+')

%% LB1null vs WT, Perimeter vs Area
h = figure;
x = 0:0.01:3.5;
scatter(stats(20).area(1).all.data,stats(20).perimeter(1).all.data,'Marker','o')
hold on;
scatter(stats(7).area(1).all.data,stats(7).perimeter(1).all.data,'Marker','+')
line(x,sqrt(x/pi).*2*pi,'Color','k')

%% LB1null vs WT, Average edge length vs Area
h = figure;
scatter(stats(20).area(1).all.data,stats(20).perimeter(1).all.data ./ stats(20).edgesPerFace(1).all.data,'Marker','o')
hold on;
scatter(stats(7).area(1).all.data,stats(7).perimeter(1).all.data ./ stats(7).edgesPerFace(1).all.data,'Marker','+')

%% LB1null vs WT, Eccentricity vs Area
h = figure;
scatter(stats(20).area(1).all.data,stats(20).eccentricity(1).all.data,'Marker','o')
hold on;
scatter(stats(7).area(1).all.data,stats(7).eccentricity(1).all.data,'Marker','+')




%% LACnull vs WT, Edges per face vs Area
h = figure;
scatter(stats(20).area(2).all.data,stats(20).edgesPerFace(2).all.data,'Marker','o')
hold on;
scatter(stats(2).area(1).all.data,stats(2).edgesPerFace(1).all.data,[],'g','Marker','+')
xlabel('Area (\mum^2)');
ylabel('Edges per Face');
legend({'WT','LAC null'});

%% LACnull vs WT, Perimeter vs Area
h = figure;
x = 0:0.01:3.5;
scatter(stats(20).area(2).all.data,stats(20).perimeter(2).all.data,'Marker','o')
hold on;
scatter(stats(2).area(1).all.data,stats(2).perimeter(1).all.data,[],'g','Marker','+')
line(x,sqrt(x/pi).*2*pi,'Color','k')
xlabel('Area (\mum^2)');
ylabel('Perimeter');
legend({'WT','LAC null','Circle'});


%% LACnull vs WT, Average edge length vs Area
h = figure;
scatter(stats(20).area(2).all.data,stats(20).perimeter(2).all.data ./ stats(20).edgesPerFace(2).all.data,'Marker','o')
hold on;
scatter(stats(2).area(1).all.data,stats(2).perimeter(1).all.data ./ stats(2).edgesPerFace(1).all.data,[],'g','Marker','+')
xlabel('Area (\mum^2)');
ylabel('Average edge length (\mum)');
legend({'WT','LAC null'});


%% LACnull vs WT, Eccentricity vs Area
h = figure;
scatter(stats(20).area(2).all.data,stats(20).eccentricity(2).all.data,'Marker','o')
hold on;
scatter(stats(2).area(1).all.data,stats(2).eccentricity(1).all.data,[],'g','Marker','+')
xlabel('Area (\mum^2)');
ylabel('Eccentricity');
legend({'WT','LAC null'});
