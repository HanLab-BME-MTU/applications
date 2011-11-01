
function res = plotMasterSlaveLifetimes(dataMaster, dataSlave, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataMaster');
ip.addRequired('dataSlave');
ip.addParamValue('ChannelNames', {'Master', 'Slave'});
ip.parse(dataMaster, dataSlave, varargin{:});


% add input checks to verify that all members are equal

dt = dataMaster(1).framerate;

cutidx = 3;
endidx = 100;

nd = numel(dataMaster);

t = cutidx*dt:dt:endidx*dt;

N = dataMaster(1).movieLength-2;
lftCorr = N./(N:-1:1);
lftCorr = lftCorr(cutidx:endidx);



for k = 1:nd
    k
    tracksMaster = loadTracks(dataMaster(k), 'Type', 'valid', 'Cutoff', cutidx);
    tracksSlave = loadTracks(dataSlave(k), 'Type', 'valid', 'Cutoff', cutidx);
    
    % generate histograms
    %niMaster = hist([tracksMaster.lifetime_s], t);
    %niMaster = niMaster / (sum(niMaster)*dt);
    
    %niSlave = hist([tracksSlave.lifetime_s], t);
    %niSlave = niSlave / (sum(niSlave)*dt);
    
    % indexes of master tracks with slave signal, and corresponding slave tracks    
    idx = assignTracksToMaster(tracksSlave, tracksMaster, 'MinOverlap', 1);    
    validMasterIndex = find(cellfun(@(c) ~isempty(c), idx));
    validSlaveIndex = idx(validMasterIndex);
    validSlaveIndex = vertcat(validSlaveIndex{:});
    
    niMaster = hist([tracksMaster(validMasterIndex).lifetime_s], t);
    niMaster = niMaster .* lftCorr;
    res.niMaster{k} = niMaster / (sum(niMaster)*dt);
    res.nMasterTracks(k) = numel(validMasterIndex);
    
    niSlave = hist([tracksSlave(validSlaveIndex).lifetime_s], t);
    niSlave = niSlave .* lftCorr;
    res.niSlave{k} = niSlave / (sum(niSlave)*dt);
    res.nSlaveTracks(k) = numel(validSlaveIndex);
    
end

% Mean histograms, SEM
M = vertcat(res.niMaster{:});
res.histMaster_Mean = mean(M,1);
res.histMaster_SEM = std(M, [], 1) / sqrt(nd);

S = vertcat(res.niSlave{:});
res.histSlave_Mean = mean(S,1);
res.histSlave_SEM = std(S, [], 1) / sqrt(nd);

res.t = t;

hueMaster = getFluorophoreHues(dataMaster(1).markers{1});
hueSlave = getFluorophoreHues(dataSlave(1).markers{1});

masterColor = hsv2rgb([hueMaster 1 0.8]);
slaveColor = hsv2rgb([hueSlave 1 0.8]);

masterFillColor = hsv2rgb([hueMaster 0.4 1]);
slaveFillColor = hsv2rgb([hueSlave 0.4 1]);

tfont = {'FontName', 'Helvetica', 'FontSize', 16};
sfont = {'FontName', 'Helvetica', 'FontSize', 20};
lfont = {'FontName', 'Helvetica', 'FontSize', 24};


figure('Position', [440 378 560 360], 'PaperPositionMode', 'auto');
hold on;
fill([t t(end:-1:1)], [res.histSlave_Mean+res.histSlave_SEM res.histSlave_Mean(end:-1:1)-res.histSlave_SEM(end:-1:1)],...
    slaveFillColor, 'EdgeColor', 'none');
hp(2) = plot(res.t, res.histSlave_Mean, '-', 'Color', slaveColor, 'LineWidth', 3);
fill([t t(end:-1:1)], [res.histMaster_Mean+res.histMaster_SEM res.histMaster_Mean(end:-1:1)-res.histMaster_SEM(end:-1:1)],...
    masterFillColor, 'EdgeColor', 'none');
hp(1) = plot(res.t, res.histMaster_Mean, '-', 'Color', masterColor, 'LineWidth', 3);
axis([0 150 0 0.05]);
set(gca, 'LineWidth', 2, 'YTick', 0:0.01:0.05, 'Layer', 'top', sfont{:});
xlabel('Lifetime (s)', lfont{:});
ylabel('Frequency', lfont{:});
hl = legend(hp, ip.Results.ChannelNames);
set(hl, 'Box', 'off');

