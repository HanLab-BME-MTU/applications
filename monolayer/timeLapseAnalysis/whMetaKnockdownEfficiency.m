function [] = whMetaKnockdownEfficiency(metaData,strLabels,outDirname)
close all;
fontsize = 24;

uniqueGenes = whGetUniqueGeneLabels(strLabels);
N = length(uniqueGenes);

strPSup = 'pSuper';
strNT = 'NT';

nSeq = 0;
allSeqStr = {};
allSeqInds = {};
allSeqKD = [];
for gene = 1 : N        
        if strcmp(uniqueGenes{gene},strPSup) || strcmp(uniqueGenes{gene},strNT)
            continue;
        end
        
        indsGene = strcmp(whToPrefix(strLabels),uniqueGenes{gene});
        
        [seqStr, seqInds] = whGetSequencesStr(strLabels,uniqueGenes{gene});
        
        for seq = 1 : length(seqStr)
            currSeqInds = seqInds{seq};            
            if sum(currSeqInds) == 0
                continue;
            end
            seqKDEff = cell2mat(metaData.KD(currSeqInds));
            seqKDEff = seqKDEff(~isnan(seqKDEff) & seqKDEff > 0);
            
            if length(seqKDEff) == 0
                continue;
            end
            
            seqKDEff = mean(seqKDEff);
           
            nSeq = nSeq + 1;
            allSeqStr{nSeq} = [uniqueGenes{gene} '_' seqStr{seq}];
            allSeqInds{nSeq} = find(currSeqInds);
            allSeqKD = [allSeqKD seqKDEff];
        end        
end

[nelements, xcenters] = hist(allSeqKD,5:10:95);
allSeqKnockdownDistribution = nelements ./ sum(nelements);

h = figure;
hold on;
bar(xcenters,allSeqKnockdownDistribution,'r');
xlabel('GEF KD efficiency (%)','FontSize',22);
ylabel('Percent','FontSize',22);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[0,100]);
set(haxes,'XTick',0:25:100);
set(haxes,'XTickLabel',0:25:100);
% set(haxes,'YLim',[0,0.13]);
% set(haxes,'YTick',0:0.05:0.1);
% set(haxes,'YTickLabel',0:0.05:0.1);
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
outFname = [outDirname 'shSeqEfficiencyDistribution.eps'];
export_fig(outFname);

h = figure;
hold on;
bar(xcenters,cumsum(allSeqKnockdownDistribution(end:-1:1)),'r');
xlabel('GEF KD efficiency (%)','FontSize',22);
ylabel('Percent','FontSize',22);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',[0,100]);
set(haxes,'XTick',0:25:100);
set(haxes,'XTickLabel',100:-25:0);
% set(haxes,'YLim',[0,1]);
% set(haxes,'YTick',0:0.2:1);
% set(haxes,'YTickLabel',0:0.2:1);
% set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
outFname = [outDirname 'shSeqEfficiencyDistributionAccumulated.eps'];
export_fig(outFname);

end