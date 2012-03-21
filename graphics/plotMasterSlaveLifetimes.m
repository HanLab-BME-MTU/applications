
% Francois Aguet (last modified: 03/20/2012)

function res = plotMasterSlaveLifetimes(dataMaster, dataSlave, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataMaster');
ip.addRequired('dataSlave');
ip.addParamValue('ChannelNames', {'Master', 'Slave'});
ip.parse(dataMaster, dataSlave, varargin{:});
chNames = ip.Results.ChannelNames;

% add input checks to verify that all members are equal


nd = numel(dataMaster);

fprintf('Master/slave plot - data sets analyzed:     ');
for k = 1:nd
    
    tracksM = loadTracks(dataMaster(k));
    tracksS = loadTracks(dataSlave(k));
    [idx unassignedSlaveIdx unassignedMasterIdx] = assignTracksToMaster(tracksS, tracksM);
    
    mappedTracksM = tracksM(cellfun(@(i) ~isempty(i), idx));
    mappedTracksS = tracksS(unique(vertcat(idx{:})));
    
    notmappedTracksM = tracksM(unassignedMasterIdx);
    notmappedTracksS = tracksS(unassignedSlaveIdx);
    
    tmp = getLifetimeHistogram(dataMaster(k), mappedTracksM);
    res.lftHistM_mapped_Ia{k} = tmp.Ia;
    
    tmp = getLifetimeHistogram(dataSlave(k), mappedTracksS);
    res.lftHistS_mapped_Ia{k} = tmp.Ia;
    
    tmp = getLifetimeHistogram(dataMaster(k), notmappedTracksM);
    res.lftHistM_notmapped_Ia{k} = tmp.Ia;
    
    tmp = getLifetimeHistogram(dataSlave(k), notmappedTracksS);
    res.lftHistS_notmapped_Ia{k} = tmp.Ia;
    
    res.pctMappedM(k) = numel(mappedTracksM)/numel(tracksM);
    res.pctMappedS(k) = numel(mappedTracksS)/numel(tracksS);
    
    % Difference in master/slave signal appearance/disappearance (class Ia tracks)
    validGaps = arrayfun(@(t) max([t.gapStatus 4]), tracksM)==4;
    idxIa_mapped = find(cellfun(@(i) ~isempty(i), idx)' & [tracksM.nSeg]==1 & [tracksM.visibility]==1 & validGaps);

    N = numel(idxIa_mapped);
    res.dStart{k} = NaN(1,N);
    res.dEnd{k} = NaN(1,N);
    for i = 1:N
        s0 = tracksM(idxIa_mapped(i)).start;
        e0 = tracksM(idxIa_mapped(i)).end;
        slaveIdx = idx{idxIa_mapped(i)};
        sSlave = min([tracksS(slaveIdx).start]);
        eSlave = max([tracksS(slaveIdx).end]);
        res.dStart{k}(i) = sSlave-s0;
        res.dEnd{k}(i) = eSlave-e0;
    end
    
    fprintf('\b\b\b\b%3d%%', round(100*k/nd));
end
fprintf('\n');

% lftHist.t = t;
% lftHist.Ia = lftHist_Ia;
% lftHist.Ib = lftHist_Ib;
% lftHist.Iab = lftHist_Iab;
% lftHist.IIa = lftHist_IIa;


% Mean histograms, SEM
M = vertcat(res.lftHistM_mapped_Ia{:});
res.meanLftHistM_mapped = mean(M,1);
% res.meanLftHistM_mapped_SEM = std(M, [], 1) / sqrt(nd);
M = vertcat(res.lftHistM_notmapped_Ia{:});
res.meanLftHistM_notmapped = mean(M,1);


S = vertcat(res.lftHistS_mapped_Ia{:});
res.meanLftHistS_mapped = mean(S,1);
% res.meanLftHistS_mapped_SEM = std(S, [], 1) / sqrt(nd);
S = vertcat(res.lftHistS_notmapped_Ia{:});
res.meanLftHistS_notmapped = mean(S,1);


% Assuming that all input data have same movieLength
res.t = tmp.t;




hueMaster = getFluorophoreHues(dataMaster(1).markers{1});
hueSlave = getFluorophoreHues(dataSlave(1).markers{1});

masterColor = hsv2rgb([hueMaster 1 0.8]);
slaveColor = hsv2rgb([hueSlave 1 0.8]);

% masterFillColor = hsv2rgb([hueMaster 0.4 1]);
% slaveFillColor = hsv2rgb([hueSlave 0.4 1]);

tfont = {'FontName', 'Helvetica', 'FontSize', 16};
sfont = {'FontName', 'Helvetica', 'FontSize', 20};
lfont = {'FontName', 'Helvetica', 'FontSize', 24};


figure('Position', [440 378 560 360], 'PaperPositionMode', 'auto');
hold on;
% fill([res.t res.t(end:-1:1)], [res.histSlave_Mean+res.histSlave_SEM res.histSlave_Mean(end:-1:1)-res.histSlave_SEM(end:-1:1)],...
%     slaveFillColor, 'EdgeColor', 'none');
hp(2) = plot(res.t, res.meanLftHistS_mapped, '.-', 'Color', slaveColor, 'LineWidth', 3, 'MarkerSize', 20);
% fill([res.t res.t(end:-1:1)], [res.histMaster_Mean+res.histMaster_SEM res.histMaster_Mean(end:-1:1)-res.histMaster_SEM(end:-1:1)],...
%     masterFillColor, 'EdgeColor', 'none');
hp(1) = plot(res.t, res.meanLftHistM_mapped, '.-', 'Color', masterColor, 'LineWidth', 3, 'MarkerSize', 20);
axis([0 140 0 0.05]);
set(gca, 'LineWidth', 2, 'YTick', 0:0.01:0.05, 'XTick', 0:20:140, 'Layer', 'top', sfont{:}, 'TickDir', 'out');
xlabel('Lifetime (s)', lfont{:});
ylabel('Frequency', lfont{:});
hl = legend(hp, chNames);
set(hl, 'Box', 'off');




figure('Position', [440 378 560 360], 'PaperPositionMode', 'auto');
hold on;
% fill([res.t res.t(end:-1:1)], [res.histSlave_Mean+res.histSlave_SEM res.histSlave_Mean(end:-1:1)-res.histSlave_SEM(end:-1:1)],...
%     slaveFillColor, 'EdgeColor', 'none');
hp(2) = plot(res.t, res.meanLftHistS_notmapped, '.-', 'Color', slaveColor, 'LineWidth', 3, 'MarkerSize', 20);
% fill([res.t res.t(end:-1:1)], [res.histMaster_Mean+res.histMaster_SEM res.histMaster_Mean(end:-1:1)-res.histMaster_SEM(end:-1:1)],...
%     masterFillColor, 'EdgeColor', 'none');
hp(1) = plot(res.t, res.meanLftHistM_notmapped, '.-', 'Color', masterColor, 'LineWidth', 3, 'MarkerSize', 20);
axis([0 140 0 0.05]);
set(gca, 'LineWidth', 2, 'YTick', 0:0.01:0.05, 'XTick', 0:20:140, 'Layer', 'top', sfont{:}, 'TickDir', 'out');
xlabel('Lifetime (s)', lfont{:});
ylabel('Frequency', lfont{:});
hl = legend(hp, chNames);
set(hl, 'Box', 'off');


dStart = [res.dStart{:}];
dStart(dStart>40) = [];
dEnd = [res.dEnd{:}];
dEnd(dEnd<-40) = [];

figure; hist(dStart, -10:40);
figure; hist(dEnd, -40:10);
