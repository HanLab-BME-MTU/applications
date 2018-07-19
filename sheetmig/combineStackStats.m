function [aveStats]=combineStackStats(groupedStats)
minLength=Inf;
for k=1:length(groupedStats)
   minLength=min(minLength,length(groupedStats(k).allStats.posVec));
end

% bring all to the same minimal length:
for k=1:length(groupedStats)
   groupedStats(k).allStats.stats =groupedStats(k).allStats.stats(1:minLength);
   groupedStats(k).allStats.posVec=groupedStats(k).allStats.posVec(1:minLength);
   groupedStats(k).allStats.numNuc=groupedStats(k).allStats.numNuc(1:minLength);
end
aveStats.posVec=groupedStats(1).allStats.posVec;

% calculate the average values:
aveStats.Ippix_MLC   =[];
aveStats.Ippix_Actin =[];
aveStats.Ippix_Nuclei=[];

aveStats.Ippix_MLC_full   =[];
aveStats.Ippix_Actin_full =[];
aveStats.Ippix_Nuclei_full=[];

aveStats.MLC_full_per_cell   =[];
aveStats.Actin_full_per_cell =[];

for k=1:length(groupedStats)
    aveStats.Ippix_MLC   =vertcat(aveStats.Ippix_MLC   ,horzcat(groupedStats(k).allStats.stats.Ippix_MLC   ));
    aveStats.Ippix_Actin =vertcat(aveStats.Ippix_Actin ,horzcat(groupedStats(k).allStats.stats.Ippix_Actin ));
    aveStats.Ippix_Nuclei=vertcat(aveStats.Ippix_Nuclei,horzcat(groupedStats(k).allStats.stats.Ippix_Nuclei));
    
    aveStats.Ippix_MLC_full   =vertcat(aveStats.Ippix_MLC_full   ,horzcat(groupedStats(k).allStats.stats.Ippix_MLC_full   ));
    aveStats.Ippix_Actin_full =vertcat(aveStats.Ippix_Actin_full ,horzcat(groupedStats(k).allStats.stats.Ippix_Actin_full ));
    aveStats.Ippix_Nuclei_full=vertcat(aveStats.Ippix_Nuclei_full,horzcat(groupedStats(k).allStats.stats.Ippix_Nuclei_full));
    
    aveStats.MLC_full_per_cell   =vertcat(aveStats.MLC_full_per_cell   ,horzcat(groupedStats(k).allStats.stats.Ippix_MLC_full  )./groupedStats(k).allStats.numNuc);
    aveStats.Actin_full_per_cell =vertcat(aveStats.Actin_full_per_cell ,horzcat(groupedStats(k).allStats.stats.Ippix_Actin_full)./groupedStats(k).allStats.numNuc);
end
aveStats.Ippix_MLC_std     =  std(aveStats.Ippix_MLC);
aveStats.Ippix_MLC_mean    = mean(aveStats.Ippix_MLC);
aveStats.Ippix_Actin_std   =  std(aveStats.Ippix_Actin);
aveStats.Ippix_Actin_mean  = mean(aveStats.Ippix_Actin);
aveStats.Ippix_Nuclei_std  =  std(aveStats.Ippix_Nuclei);
aveStats.Ippix_Nuclei_mean = mean(aveStats.Ippix_Nuclei);

aveStats.Ippix_MLC_full_std     =  std(aveStats.Ippix_MLC_full);
aveStats.Ippix_MLC_full_mean    = mean(aveStats.Ippix_MLC_full);
aveStats.Ippix_Actin_full_std   =  std(aveStats.Ippix_Actin_full);
aveStats.Ippix_Actin_full_mean  = mean(aveStats.Ippix_Actin_full);
aveStats.Ippix_Nuclei_full_std  =  std(aveStats.Ippix_Nuclei_full);
aveStats.Ippix_Nuclei_full_mean = mean(aveStats.Ippix_Nuclei_full);

aveStats.MLC_full_per_cell_std    =  std(aveStats.MLC_full_per_cell);
aveStats.MLC_full_per_cell_mean   = mean(aveStats.MLC_full_per_cell);
aveStats.Actin_full_per_cell_std  =  std(aveStats.Actin_full_per_cell);
aveStats.Actin_full_per_cell_mean = mean(aveStats.Actin_full_per_cell);


figure()
errorbar(aveStats.posVec,aveStats.Ippix_MLC_mean,aveStats.Ippix_MLC_std);
title('average MLC-Intensity per pixel');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'figures3kPa/MLC.fig');

figure()
errorbar(aveStats.posVec,aveStats.Ippix_Actin_mean,aveStats.Ippix_Actin_std);
title('average Actin-Intensity per pixel');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'figures3kPa/Actin.fig');

figure()
errorbar(aveStats.posVec,aveStats.Ippix_Nuclei_mean,aveStats.Ippix_Nuclei_std);
title('average Nuclei-Intensity per pixel');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'figures3kPa/Nuclei.fig');

% MLC-Intensity normalized:

figure()
% the std is:
errorbar(aveStats.posVec,aveStats.Ippix_MLC_mean./aveStats.Ippix_Actin_mean,aveStats.Ippix_MLC_std./aveStats.Ippix_Actin_mean);
title('tot. MLC-Int. normalized by tot Actin-Int.');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'figures3kPa/MLCOverActin.fig');

figure()
errorbar(aveStats.posVec,aveStats.Ippix_MLC_mean./aveStats.Ippix_Nuclei_mean,aveStats.Ippix_MLC_std./aveStats.Ippix_Nuclei_mean);
title('tot. MLC-Int. normalized by tot Nuclei-Int.');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'figures3kPa/MLCOverNuclei.fig');

figure()
errorbar(aveStats.posVec,aveStats.Ippix_Actin_mean./aveStats.Ippix_Nuclei_mean,aveStats.Ippix_Actin_std./aveStats.Ippix_Nuclei_mean);
title('tot. Actin-Int. normalized by tot Nuclei-Int.');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'figures3kPa/ActinOverNuclei.fig');


% MLC-Intensity normalized by number of cells:
%return;

figure()
%errorbar(posVec,horzcat(stats.Ippix_MLC_full)./numNuc,facSEMtoSEM95*horzcat(stats.Ippix_MLC_std)./numNuc./sqrt(numNuc));
errorbar(aveStats.posVec,aveStats.MLC_full_per_cell_mean,aveStats.MLC_full_per_cell_std);
title('MLC-Int normalized by Nuclei Number');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'figures3kPa/MLCOverNucleiNumber.fig');

figure()
%errorbar(posVec,horzcat(stats.Ippix_Actin_full)./numNuc,facSEMtoSEM95*horzcat(stats.Ippix_Actin_std)./numNuc./sqrt(numNuc));
errorbar(aveStats.posVec,aveStats.Actin_full_per_cell_mean,aveStats.Actin_full_per_cell_std);
title('Actin-Int. normalized by Nuclei Number');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
box on
set(gca,'FontSize',18)
saveas(gcf,'figures3kPa/ActinOverNucleiNumber.fig');


return;
figure()
errorbar(aveStats3kpa.posVec,aveStats3kpa.Ippix_MLC_mean,aveStats3kpa.Ippix_MLC_std,'-sb');
hold on;
errorbar(aveStats65kpa.posVec,aveStats65kpa.Ippix_MLC_mean,aveStats65kpa.Ippix_MLC_std,'-sr');
hold off;
title('average MLC-Intensity per pixel');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
xlim([50 750])
box on
set(gca,'FontSize',18)
saveas(gcf,'figuresComp/MLC.fig');
saveas(gcf,'figuresComp/MLC.eps','psc2');


% compare soft and stiff:
figure()
% the std is:
errorbar(aveStats3kpa.posVec,aveStats3kpa.Ippix_MLC_mean./aveStats3kpa.Ippix_Actin_mean,aveStats3kpa.Ippix_MLC_std./aveStats3kpa.Ippix_Actin_mean,'-sb');
hold on;
errorbar(aveStats65kpa.posVec,aveStats65kpa.Ippix_MLC_mean./aveStats65kpa.Ippix_Actin_mean,aveStats65kpa.Ippix_MLC_std./aveStats65kpa.Ippix_Actin_mean,'-sr');
hold off
title('tot. MLC-Int. normalized by tot Actin-Int.');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
xlim([50 750])
box on
set(gca,'FontSize',18)
saveas(gcf,'figuresComp/MLCOverActin.fig');
saveas(gcf,'figuresComp/MLCOverActin.eps','psc2');


figure()
%errorbar(posVec,horzcat(stats.Ippix_MLC_full)./numNuc,facSEMtoSEM95*horzcat(stats.Ippix_MLC_std)./numNuc./sqrt(numNuc));
errorbar(aveStats3kpa.posVec,aveStats3kpa.MLC_full_per_cell_mean,aveStats3kpa.MLC_full_per_cell_std,'-sb');
hold on;
errorbar(aveStats65kpa.posVec,aveStats65kpa.MLC_full_per_cell_mean,aveStats65kpa.MLC_full_per_cell_std,'-sr');
hold off
title('MLC-Int normalized by Nuclei Number');
xlabel('position relative to sheet edge [um]');
ylabel('I [a.u.]');
xlim([50 750])
box on
set(gca,'FontSize',18)
saveas(gcf,'figuresComp/MLCOverNucleiNumber.fig');
saveas(gcf,'figuresComp/MLCOverNucleiNumber.eps','psc2');
