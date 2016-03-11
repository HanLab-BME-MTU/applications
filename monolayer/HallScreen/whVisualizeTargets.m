
% targetGenesStr - cell array of strings for genes of interest
function [] = whVisualizeTargets(geneDayDiff,allDiffMeanGeneToMeanControl,targetGenesStr,mainDirname,propertyStr)



[coeffMeanGeneMeanControl,scoreMeanGeneMeanControl,latentMeanGeneMeanControl] = pca(allDiffMeanGeneToMeanControl');
accVarianceMeanGeneMeanControl = cumsum(latentMeanGeneMeanControl)./sum(latentMeanGeneMeanControl);

plotExperiments(geneDayDiff,targetGenesStr,mainDirname,propertyStr,coeffMeanGeneMeanControl,scoreMeanGeneMeanControl,latentMeanGeneMeanControl,2);
plotExperiments(geneDayDiff,targetGenesStr,mainDirname,propertyStr,coeffMeanGeneMeanControl,scoreMeanGeneMeanControl,latentMeanGeneMeanControl,3);
end

function plotExperiments(geneDayDiff,targetGenesStr,mainDirname,propertyStr,coeffMeanGeneMeanControl,scoreMeanGeneMeanControl,latentMeanGeneMeanControl,PCn)

nTargets = length(targetGenesStr);
nConditions = nTargets + 3; % 0% KD, CDC42/RAC1/beta-PIX, > 50% KD (rest)

cmap = colormap(hsv(nConditions));
fontsize = 24;
markerSize = 7;
LineWidth = 3;


[negCtrlInds,restInds,targetsInds,posCntrl] = whGetTargetInds(geneDayDiff,targetGenesStr);

nGeneDay = length(geneDayDiff);


% negCtrlInds = [];
% restInds = [];
% targetsInds = cell(1,nTargets);
% posCntrl = [];
% 
% 
% for i = 1 : nTargets
%     targetsInds{i} = [];
% end
% 
% for iGeneDay = 1 : nGeneDay  
%     geneStr = geneDayDiff{iGeneDay}.geneStr;
%     KD = geneDayDiff{iGeneDay}.KD;
%     nameStr = nan;
%     if KD == 0
%        negCtrlInds = [negCtrlInds iGeneDay];
%        continue;
%     else if strcmp(geneStr,'CDC42') || strcmp(geneStr,'RAC1') || strcmp(geneStr,'beta-PIX')
%             posCntrl = [posCntrl iGeneDay];
%             continue;
%         else
%             check = false;
%             for t = 1 : nTargets
%                 if strcmp(geneStr,targetGenesStr{t})
%                     targetsInds{t} = [targetsInds{t} iGeneDay];
%                     check = true;
%                     break;
%                 end
%             end            
%             % > 50%
%             if ~check
%                 restInds = [restInds iGeneDay];
%             end
%         end
%     end    
% end



h = figure;
xlabel('PC 1','FontSize',fontsize);
ylabel(['PC ' num2str(PCn)],'FontSize',fontsize);
hold on;

plot(scoreMeanGeneMeanControl(negCtrlInds,1),scoreMeanGeneMeanControl(negCtrlInds,PCn),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','0% KD');   
plot(scoreMeanGeneMeanControl(restInds,1),scoreMeanGeneMeanControl(restInds,PCn),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','> 50% KD');   
for t = 1 : nTargets
% for t = 4 : 4
plot(scoreMeanGeneMeanControl(targetsInds{t},1),scoreMeanGeneMeanControl(targetsInds{t},PCn),'o','MarkerEdgeColor',cmap(t+2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName',targetGenesStr{t});
%     plot(scoreMeanGeneMeanControl(targetsInds{t},1),scoreMeanGeneMeanControl(targetsInds{t},PCn),'o','MarkerEdgeColor',cmap(t+2,:),'LineWidth',LineWidth,'MarkerSize',markerSize+3,'DisplayName',targetGenesStr{t});
end
plot(scoreMeanGeneMeanControl(posCntrl.inds,1),scoreMeanGeneMeanControl(posCntrl.inds,PCn),'o','MarkerEdgeColor',cmap(nConditions,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','CDC42/RAC1/beta-PIX');

haxes = get(h,'CurrentAxes');
if strcmp(propertyStr,'Speed')
    assert(max(abs(scoreMeanGeneMeanControl(:,1))) < 100);
    set(haxes,'XLim',[-100,100]);
    set(haxes,'XTick',-80:40:80);
    set(haxes,'XTickLabel',-80:40:80);
    if PCn == 2
        assert(max(abs(scoreMeanGeneMeanControl(:,2))) < 60);
        set(haxes,'YLim',[-60,60]);
        set(haxes,'YTick',-50:25:50);
        set(haxes,'YTickLabel',-50:25:50);
    else if PCn == 3
            assert(max(abs(scoreMeanGeneMeanControl(:,3))) < 23);
            set(haxes,'YLim',[-23,23]);
            set(haxes,'YTick',-20:10:20);
            set(haxes,'YTickLabel',-20:10:20);
        else
            error(['PC ' nPC ' > 3']);
        end
    end
else if strcmp(propertyStr,'Directionality')
        assert(max(abs(scoreMeanGeneMeanControl(:,1))) < 9);
        set(haxes,'XLim',[-9,9]);
        set(haxes,'XTick',-8:4:8);
        set(haxes,'XTickLabel',-8:4:8);
        if PCn == 2
            assert(max(abs(scoreMeanGeneMeanControl(:,2))) < 5);
            set(haxes,'YLim',[-5,5]);
            set(haxes,'YTick',-5:2.5:5);
            set(haxes,'YTickLabel',-5:2.5:5);
        else if PCn == 3
                assert(max(abs(scoreMeanGeneMeanControl(:,3))) < 3.2);
                set(haxes,'YLim',[-3.2,3.2]);
                set(haxes,'YTick',-3:1.5:3);
                set(haxes,'YTickLabel',-3:1.5:3);
            else
                error(['PC ' nPC ' > 3']);
            end
        end
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(scoreMeanGeneMeanControl(:,1))) < 1.2);
            set(haxes,'XLim',[-1.2,1.2]);
            set(haxes,'XTick',-1:0.5:1);
            set(haxes,'XTickLabel',-1:0.5:1);
            if PCn == 2
                assert(max(abs(scoreMeanGeneMeanControl(:,2))) < 0.8);
                set(haxes,'YLim',[-0.8,0.8]);
                set(haxes,'YTick',-0.8:0.4:0.8);
                set(haxes,'YTickLabel',-0.8:0.4:0.8);
            else if PCn == 3
                    assert(max(abs(scoreMeanGeneMeanControl(:,3))) < 0.42);
                    set(haxes,'YLim',[-0.42,0.42]);
                    set(haxes,'YTick',-0.4:0.2:0.4);
                    set(haxes,'YTickLabel',-0.4:0.2:0.4);
                else
                    error(['PC ' nPC ' > 3']);
                end
            end
            
        end
    end
end

set(haxes,'FontSize',fontsize);

legend('show');

plot(0,0,'*k','LineWidth',4,'MarkerSize',20); 

set(h,'Color','none');

controlGeneVarFname = [mainDirname 'targets/' propertyStr '_PC1PC' num2str(PCn) '_legend.eps'];
export_fig(controlGeneVarFname);

legend off;
controlGeneVarFname = [mainDirname 'targets/' propertyStr '_PC1PC' num2str(PCn) '.eps'];
export_fig(controlGeneVarFname);

hold off;

%% Positive control plot
h = figure;
xlabel('PC 1','FontSize',fontsize);
ylabel(['PC ' num2str(PCn)],'FontSize',fontsize);
hold on;

plot(scoreMeanGeneMeanControl(negCtrlInds,1),scoreMeanGeneMeanControl(negCtrlInds,PCn),'o','MarkerEdgeColor',cmap(1,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','0% KD');   
plot(scoreMeanGeneMeanControl(restInds,1),scoreMeanGeneMeanControl(restInds,PCn),'o','MarkerEdgeColor',cmap(2,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName','> 50% KD');   
for t = 1 : 3
    plot(scoreMeanGeneMeanControl(posCntrl.targetInds{t},1),scoreMeanGeneMeanControl(posCntrl.targetInds{t},PCn),'o','MarkerEdgeColor',cmap(2+t,:),'LineWidth',LineWidth,'MarkerSize',markerSize,'DisplayName',posCntrl.targetStrs{t});                
end

haxes = get(h,'CurrentAxes');
if strcmp(propertyStr,'Speed')
    assert(max(abs(scoreMeanGeneMeanControl(:,1))) < 100);
    set(haxes,'XLim',[-100,100]);
    set(haxes,'XTick',-80:40:80);
    set(haxes,'XTickLabel',-80:40:80);
    if PCn == 2
        assert(max(abs(scoreMeanGeneMeanControl(:,2))) < 60);
        set(haxes,'YLim',[-60,60]);
        set(haxes,'YTick',-50:25:50);
        set(haxes,'YTickLabel',-50:25:50);
    else if PCn == 3
            assert(max(abs(scoreMeanGeneMeanControl(:,3))) < 23);
            set(haxes,'YLim',[-23,23]);
            set(haxes,'YTick',-20:10:20);
            set(haxes,'YTickLabel',-20:10:20);
        else
            error(['PC ' nPC ' > 3']);
        end
    end
else if strcmp(propertyStr,'Directionality')
        assert(max(abs(scoreMeanGeneMeanControl(:,1))) < 9);
        set(haxes,'XLim',[-9,9]);
        set(haxes,'XTick',-8:4:8);
        set(haxes,'XTickLabel',-8:4:8);
        if PCn == 2
            assert(max(abs(scoreMeanGeneMeanControl(:,2))) < 5);
            set(haxes,'YLim',[-5,5]);
            set(haxes,'YTick',-5:2.5:5);
            set(haxes,'YTickLabel',-5:2.5:5);
        else if PCn == 3
                assert(max(abs(scoreMeanGeneMeanControl(:,3))) < 3.2);
                set(haxes,'YLim',[-3.2,3.2]);
                set(haxes,'YTick',-3:1.5:3);
                set(haxes,'YTickLabel',-3:1.5:3);
            else
                error(['PC ' nPC ' > 3']);
            end
        end
    else if strcmp(propertyStr,'Coordination')
            assert(max(abs(scoreMeanGeneMeanControl(:,1))) < 1.2);
            set(haxes,'XLim',[-1.2,1.2]);
            set(haxes,'XTick',-1:0.5:1);
            set(haxes,'XTickLabel',-1:0.5:1);
            if PCn == 2
                assert(max(abs(scoreMeanGeneMeanControl(:,2))) < 0.8);
                set(haxes,'YLim',[-0.8,0.8]);
                set(haxes,'YTick',-0.8:0.4:0.8);
                set(haxes,'YTickLabel',-0.8:0.4:0.8);
            else if PCn == 3
                    assert(max(abs(scoreMeanGeneMeanControl(:,3))) < 0.42);
                    set(haxes,'YLim',[-0.42,0.42]);
                    set(haxes,'YTick',-0.4:0.2:0.4);
                    set(haxes,'YTickLabel',-0.4:0.2:0.4);
                else
                    error(['PC ' nPC ' > 3']);
                end
            end
            
        end
    end
end

set(haxes,'FontSize',fontsize);

legend('show');

plot(0,0,'*k','LineWidth',4,'MarkerSize',20); 

set(h,'Color','none');

controlGeneVarFname = [mainDirname 'targets/' propertyStr '_PC1PC' num2str(PCn) '_posCtrl_legend.eps'];
export_fig(controlGeneVarFname);

legend off;
controlGeneVarFname = [mainDirname 'targets/' propertyStr '_PC1PC' num2str(PCn) '_posCtrl.eps'];
export_fig(controlGeneVarFname);

hold off;


%% Plot PCs
pcPlot(coeffMeanGeneMeanControl,latentMeanGeneMeanControl,propertyStr,[mainDirname 'targets/']);
close all;
end