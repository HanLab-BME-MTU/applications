%% LA
% Using LA / LB1

% quantiles = 100 - [(sqrt(2).^-(0:50))*100 0];
% quantiles = [10:10:90];

h = figure;
subplot(2,2,1);
qh(1) = lamins.plot.qqplot_enhanced(stats(21).area(1).all.data,stats(29).area.all.data,quantiles);
ylabel('mEmerald-Lamin A, Face Area (\mum^2)');
xlabel('Lamin A Immunofluorescence, Face Area(\mum^2)');
%% LB1
% Using LA / LB1
% figure;
subplot(2,2,2);
qh(2) = lamins.plot.qqplot_enhanced(stats(21).area(2).all.data,stats(30).area.all.data,quantiles);
ylabel('mEmerald-Lamin B1, Face Area (\mum^2)');
xlabel('Lamin B1 Immunofluorescence, Face Area (\mum^2)');
%% LB2
% Using LC / LB2
% figure;
subplot(2,2,4);
qh(3) = lamins.plot.qqplot_enhanced(stats(28).area(2).all.data,stats(31).area.all.data,quantiles);
ylabel('mEmerald-Lamin B2, Face Area (\mum^2)');
xlabel('Lamin B2 Immunofluorescence, Face Area (\mum^2)');
%% LC
% Using LC / LB2
% figure;
subplot(2,2,3);
qh(4) = lamins.plot.qqplot_enhanced(stats(28).area(1).all.data,stats(32).area.all.data,quantiles);
ylabel('mEmerald-Lamin C, Face Area (\mum^2)');
xlabel('Lamin C Immunofluorescence, Face area (\mum^2)');
set(h,'Position',[680 355 992 742]);
set(h,'Position',[680 355 992 742]);
set(h,'PaperUnits','normalized');
% set(h,'PaperSize',[11 8.5]);
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition',[0 0 1 1])
print(h,'figS2','-dpdf');

%%% Inset Figure S2

%% LA
% Using LA / LB1

quantiles = [10:10:90];

h = figure;
subplot(2,2,1);
qh(1) = lamins.plot.qqplot_enhanced(stats(21).area(1).all.data,stats(29).area.all.data,quantiles);
ylim([0 1]);
xlim([0 1]);
ylabel('mEmerald-Lamin A, Face Area (um^2)');
xlabel('Lamin A Immunofluorescence, Face Area(um^2)');
%% LB1
% Using LA / LB1
% figure;
subplot(2,2,2);
qh(2) = lamins.plot.qqplot_enhanced(stats(21).area(2).all.data,stats(30).area.all.data,quantiles);
ylabel('mEmerald-Lamin B1, Face Area (um^2)');
xlabel('Lamin B1 Immunofluorescence, Face Area (um^2)');
ylim([0 1]);
xlim([0 1]);
%% LB2
% Using LC / LB2
% figure;
subplot(2,2,4);
qh(3) = lamins.plot.qqplot_enhanced(stats(28).area(2).all.data,stats(31).area.all.data,quantiles);
ylabel('mEmerald-Lamin B2, Face Area (um^2)');
xlabel('Lamin B2 Immunofluorescence, Face Area (um^2)');
ylim([0 1]);
xlim([0 1]);

%% LC
% Using LC / LB2
% figure;
subplot(2,2,3);
qh(4) = lamins.plot.qqplot_enhanced(stats(28).area(1).all.data,stats(32).area.all.data,quantiles);
ylabel('mEmerald-Lamin C, Face Area (um^2)');
xlabel('Lamin C Immunofluorescence, Face area (um^2)');
ylim([0 1]);
xlim([0 1]);

set(h,'Position',[680 355 992 742]);
set(h,'PaperUnits','normalized');
% set(h,'PaperSize',[11 8.5]);
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition',[0 0 1 1]);
print(h,'fig2','-dpdf');

%% LA
% Using LA / LB1

quantiles = [0:10:100];

h = figure;
subplot(2,2,1);
X = stats(21).area(1).all.data;
Y = stats(29).area.all.data;
qh(1) = lamins.plot.qqplot_enhanced(X,Y,quantiles);
Xp = prctile(X,100 - 0.5.^(1:10));
Yp = prctile(Y,100 - 0.5.^(1:10));
hold on;
scatter(Xp,Yp,'g','Marker','+');
% ylim([0 1]);
% xlim([0 1]);
ylabel('mEmerald-Lamin A, Face Area (um^2)');
xlabel('Lamin A Immunofluorescence, Face Area(um^2)');
%% LB1
% Using LA / LB1
% figure;
subplot(2,2,2);
X = stats(21).area(2).all.data;
Y = stats(30).area.all.data;
qh(2) = lamins.plot.qqplot_enhanced(X,Y,quantiles);
Xp = prctile(X,100 - 0.5.^(1:10));
Yp = prctile(Y,100 - 0.5.^(1:10));
hold on;
scatter(Xp,Yp,'g','Marker','+');
ylabel('mEmerald-Lamin B1, Face Area (um^2)');
xlabel('Lamin B1 Immunofluorescence, Face Area (um^2)');
% ylim([0 1]);
% xlim([0 1]);
%% LB2
% Using LC / LB2
% figure;
subplot(2,2,4);
qh(3) = lamins.plot.qqplot_enhanced(stats(28).area(2).all.data,stats(31).area.all.data,quantiles);
X = stats(28).area(2).all.data;
Y = stats(31).area.all.data;
Xp = prctile(X,100 - 0.5.^(1:10));
Yp = prctile(Y,100 - 0.5.^(1:10));
hold on;
scatter(Xp,Yp,'g','Marker','+');
ylabel('mEmerald-Lamin B2, Face Area (um^2)');
xlabel('Lamin B2 Immunofluorescence, Face Area (um^2)');
% ylim([0 1]);
% xlim([0 1]);

%% LC
% Using LC / LB2
% figure;
subplot(2,2,3);
qh(4) = lamins.plot.qqplot_enhanced(stats(28).area(1).all.data,stats(32).area.all.data,quantiles);
X = stats(28).area(1).all.data;
Y = stats(32).area.all.data;
Xp = prctile(X,100 - 0.5.^(1:10));
Yp = prctile(Y,100 - 0.5.^(1:10));
hold on;
scatter(Xp,Yp,'g','Marker','+');
ylabel('mEmerald-Lamin C, Face Area (um^2)');
xlabel('Lamin C Immunofluorescence, Face area (um^2)');
% ylim([0 1]);
% xlim([0 1]);

set(h,'Position',[680 355 992 742]);
set(h,'PaperUnits','normalized');
% set(h,'PaperSize',[11 8.5]);
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition',[0 0 1 1]);
print(h,'fig2KJ','-dpdf');

