%% Record relative change of LBP over time for every single cell in an experiment
% output at the corresponding analysis/Cells directory 
% NOT SUPPORTING BCK,FWD!

% Assaf Zaritsky, June. 2016

function [] = pcCellsLBP()

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/'));

close all;

nScales = 4;
scales = 1.0./2.^((1:nScales)-1);

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';

metaDataFname = [analysisDirname 'MetaData/Experiments20151023.mat'];

load(metaDataFname);%metaData

for iScale = 1 : nScales    
    visCellsLbpPlasticity(metaData,analysisDirname,iScale);        
end

end

%%
function [] = visCellsLbpPlasticity(metaData,analysisDirname,iScale)
lbpDirname = [analysisDirname 'Cells/LBP'];

fovDname = [lbpDirname filesep 'FOV' filesep]; if ~exist(fovDname,'dir') unix(sprintf('mkdir %s',fovDname)); end %#ok<SEPEX>
% bckDname = [lbpDirname filesep 'BCK' filesep]; if ~exist(bckDname,'dir') unix(sprintf('mkdir %s',bckDname)); end %#ok<SEPEX>
% fwdDname = [lbpDirname filesep 'FWD' filesep]; if ~exist(fwdDname,'dir') unix(sprintf('mkdir %s',fwdDname)); end %#ok<SEPEX>

curFovDname = [fovDname num2str(iScale) filesep]; if ~exist(curFovDname,'dir') unix(sprintf('mkdir %s',curFovDname)); end %#ok<SEPEX>
% curBckDname = [fovDname num2str(iScale) filesep]; if ~exist(curBckDname,'dir') unix(sprintf('mkdir %s',curBckDname)); end %#ok<SEPEX>
% curFwdDname = [fovDname num2str(iScale) filesep]; if ~exist(curFwdDname,'dir') unix(sprintf('mkdir %s',curFwdDname)); end %#ok<SEPEX>


for iexp = 1 : metaData.experiments.N
    curFname = metaData.experiments.fnames{iexp};
    curDate = curFname(1:6);
    
    outFovDname = [curFovDname filesep curFname filesep]; if ~exist(outFovDname,'dir') unix(sprintf('mkdir %s',outFovDname)); end %#ok<SEPEX>
    %     outBckDname = [curBckDname filesep curFname]; if ~exist(outBckDname,'dir') unix(sprintf('mkdir %s',outBckDname)); end %#ok<SEPEX>
    %     outFwdDname = [curFwdDname filesep curFname]; if ~exist(outFwdDname,'dir') unix(sprintf('mkdir %s',outFwdDname)); end %#ok<SEPEX>

    % cell lbp plasticity for the whole experiment
    lbpPlasticity = cell(1,2);
    lbpPlasticity{1}.cells = {};
    lbpPlasticity{1}.well = [];
    lbpPlasticity{2}.cells = {};
    lbpPlasticity{2}.well = [];
    
    for in = 1 : 2
        %         curAll = length(out.allCells.fov.cellLbp) + 1; % # of the current experiment / well
        
        if in == 1
            tasksItr = 1 : metaData.experiments.n1{iexp};
            curSource = metaData.experiments.source1{iexp};            
        else
            tasksItr = (metaData.experiments.n1{iexp} + 1) : (metaData.experiments.n1{iexp} + metaData.experiments.n2{iexp});
            curSource = metaData.experiments.source1{iexp};
        end
                        
        for itask = tasksItr
            
            % in exclude list
            if ismember(itask,metaData.experiments.exclude{iexp});
                continue;
            end
            
            lbpFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
                curFname '_s' sprintf('%02d',itask) filesep...
                'tracking' filesep 'lbpData.mat'];
            
            cellIdTxyFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
                curFname '_s' sprintf('%02d',itask) filesep...
                'tracking' filesep 'cellIdTYX.mat'];
            
            if ~exist(lbpFname,'file')
                lbpFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
                    curFname '_s' sprintf('%d',itask) filesep...
                    'tracking' filesep 'lbpData.mat'];
                
                cellIdTxyFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
                    curFname '_s' sprintf('%d',itask) filesep...
                    'tracking' filesep 'cellIdTYX.mat'];
            end
            
            if ~exist(lbpFname,'file')
                warning(['LBP file' lbpFname ' (either 1 or 01) does not exist']);
                continue;
                error(['LBP file' lbpFname ' (either 1 or 01) does not exist']);
            end
            
            load(lbpFname); % lbpData
            load(cellIdTxyFname); % cellTXY
            
            nCurCells = length(lbpData.fov);
            lbpPlasticityWell = cell(1,nCurCells);            
                        
            for icell = 1 : nCurCells
                ts = cellTYX{icell}.ts;
                xs = cellTYX{icell}.xs;
                ys = cellTYX{icell}.ys;
                lbp = lbpData.fov{icell}.pyramidLBP{iScale}.lbp;% 47 x 10
                plasticLbp = lbp - repmat(lbp(1,:),size(lbp,1),1); % TODO: A - repmat(v,size(A,2),1)
                plasticDist = sum(abs(plasticLbp),2);
                deltaLbp = abs(plasticDist(1:end-1)-plasticDist(2:end));
                
                lbpPlasticityWell{icell}.ts = ts;
                lbpPlasticityWell{icell}.xs = xs;
                lbpPlasticityWell{icell}.ys = ys;
                lbpPlasticityWell{icell}.lbp = lbp;
                lbpPlasticityWell{icell}.plasticLbp = plasticLbp;
                lbpPlasticityWell{icell}.plasticDist = plasticDist;
                lbpPlasticityWell{icell}.deltaLbp = deltaLbp;
            end
            wellDname = [outFovDname filesep 's' sprintf('%02d',itask) filesep];
            
            plotLbpPlasticity(lbpPlasticityWell,wellDname,'',true); % last two arguments are prefix = '' and printAllCells = true
            
            nAccCells = length(lbpPlasticity{in}.cells);
            iAccCell = icell + nAccCells; % cell # in the accumulated data structure
            
            lbpPlasticity{in}.cells = [lbpPlasticity{in}.cells,lbpPlasticityWell];                
            lbpPlasticity{in}.well = [lbpPlasticity{in}.well, itask*ones(1,length(lbpPlasticityWell))];
        end
        expPrefix = sprintf('%s_%d_lbpPlasticity',curFname,in);
        % make the plots in the main directory: outFovDname
        if ~isempty(lbpPlasticity{in}.cells)
            plotLbpPlasticity(lbpPlasticity{in}.cells,curFovDname,expPrefix,false); % last two arguments are prefix = '' and printAllCells = false
        end
    end
end



end

%% 
% printAllCells - only if 
function []  = plotLbpPlasticity(cells,outDname,outPrefix,printAllCells)
close all;

nStdTH = 4; % number of STDs above mean to be considered as an event
nPrctileTH = 0.01; % % of cells that need to be above the threshold at a given time point to be considered as a QC "event"
timeWin = 15; % minutes

if ~exist(outDname,'dir') 
    unix(sprintf('mkdir %s',outDname)); 
end

nCells = length(cells);

if printAllCells
    for icell = 1 : nCells
        outCellFname = [outDname filesep outPrefix num2str(icell) '.eps'];
        if ~exist(outCellFname,'file')
            curCell = cells{icell};
            curCellMaxDelta = max(curCell.deltaLbp);
            curCellStdDelta = std(curCell.deltaLbp);
            curCellMaxDist = max(curCell.plasticDist);
            h = figure('Visible','off');
            title(sprintf('max diff = %.3f, std = %.3f (max dist = %.3f)',curCellMaxDelta,curCellStdDelta,curCellMaxDist));
            hold on;
            plot(curCell.ts,curCell.plasticDist,'--k');
            plot(curCell.ts(1),curCell.plasticDist(1),'or','MarkerSize',12,'LineWidth',3);
            set(h,'Color','w');
            %         caxis([0,1]);
            hold off;
            export_fig(outCellFname);
            close all;
        end
    end
end


outAllCellsFname = [outDname filesep outPrefix '_stats.eps'];

if ~exist(outAllCellsFname,'file')
    
    allCellsDeltaLbp = [];
    allCellsDeltaLbpStd = [];
    maxT = 0;
    
    % Plot all cells on the same plot, calculate statistics + save data
    h = figure;%('Visible','off');
    hold on;
    for icell = 1 : nCells
        curCell = cells{icell};
        plot(curCell.ts(2:end),curCell.deltaLbp,'--k');
        plot(curCell.ts(1),0,'og','MarkerSize',6,'LineWidth',2);
        
        allCellsDeltaLbp = [allCellsDeltaLbp, curCell.deltaLbp'];
        allCellsDeltaLbpStd = [allCellsDeltaLbpStd, std(curCell.deltaLbp)];
        
        maxT = max(maxT,curCell.ts(end));
    end
    
    stdAllCells = std(allCellsDeltaLbp);
    stdAllCellsTH = mean(allCellsDeltaLbp) + nStdTH*stdAllCells;
    
    ylim([0,1.1*stdAllCellsTH]);
    
    plot([0,maxT],[stdAllCellsTH,stdAllCellsTH],'b--');
    
    timeQC = zeros(1,2000);
    for icell = 1 : nCells
        curCell = cells{icell};
        inds = curCell.ts(1) + find(curCell.deltaLbp > stdAllCellsTH);
        timeQC(inds) = timeQC(inds) + 1;        
    end
    
    timeQCAcc = accTimeQC(timeQC,timeWin);
    timeQCNorm = timeQCAcc ./ nCells;
    timeEvents = find(timeQCNorm > nPrctileTH);
    if ~isempty(timeEvents)
        plot(timeEvents,stdAllCellsTH*ones(1,length(timeEvents)),'or','MarkerSize',6,'LineWidth',2,'MarkerFaceColor','r');
    end
    
    set(h,'Color','w');
    hold off;
    export_fig(outAllCellsFname);
    save([outDname filesep outPrefix '_stats.mat'],'allCellsDeltaLbp','allCellsDeltaLbpStd','nCells','stdAllCells','stdAllCellsTH','timeQC','timeQCNorm','timeEvents','maxT');
end


end

function timeQCAcc = accTimeQC(timeQC,timeWin)
halfTimeWin = floor(timeWin/2.0);
timeQC = [zeros(1,halfTimeWin) timeQC zeros(1,halfTimeWin)];

timeQCAcc = zeros(size(timeQC));

for i = (halfTimeWin + 1) : (length(timeQC) - halfTimeWin)
    timeQCAcc(i) = sum(timeQC((i-halfTimeWin):(i+halfTimeWin)));
end

timeQCAcc = timeQCAcc((halfTimeWin + 1) : (length(timeQC) - halfTimeWin));

end