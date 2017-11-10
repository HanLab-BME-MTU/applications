% scriptCompareCCWSvsSPC
% scriptCompareCCWSvsSPC compares performance of tracking algorithm in CCWS
% and SPC for one-time-pass by making synthetic bead images
%% Shift vs. Err
p=0;
nShifts = 3/0.05+1;
errFCTR = cell(nShifts,1);
errCDWS = cell(nShifts,1);
errCCWS = cell(nShifts,1);
errSPC = cell(nShifts,1);
for shift=0:0.05:3
    p=p+1;
    errFCTR{p} = testTrackStackFlow(shift,0,17,'FCTR',5);
    errCDWS{p} = testTrackStackFlow(shift,0,17,'CDWS',5);
    errCCWS{p} = testTrackStackFlow(shift,0,17,'CCWS',5);
    errSPC{p} = testTrackStackFlow(shift,0,17,'SPC',5);
end
%% Shift vs. Err
p=0;
nShifts = 3/0.05+1;
errCDWS = cell(nShifts,1);
for shift=0:0.05:3
    p=p+1;
    errCDWS{p} = testTrackStackFlow(shift,0,17,'CDWS',5);
    errCCWS{p} = testTrackStackFlow(shift,0,17,'CCWS',5);
end
%% Image making for imageJ-based CCWS
for shift=0:0.05:3
    p=p+1;
    testTrackStackFlow(shift,0,17,'CCWS',5,true);
end

%% Converting cell array to matrix
%max number determination
maxLenFCTR=0;
maxLenCDWS=0;
maxLenCCWS=0;
maxLenSPC=0;
for ii=1:nShifts
    maxLenFCTR = max(maxLenFCTR,length(errFCTR{ii}));
    maxLenCDWS = max(maxLenCDWS,length(errCDWS{ii}));
    maxLenCCWS = max(maxLenCCWS,length(errCCWS{ii}));
    maxLenSPC = max(maxLenSPC,length(errSPC{ii}));
end
errFCTRmat = NaN(maxLenFCTR,nShifts);
errCDWSmat = NaN(maxLenCDWS,nShifts);
errCCWSmat = NaN(maxLenCCWS,nShifts);
errSPCmat = NaN(maxLenSPC,nShifts);
for ii=1:nShifts
    curLenFCTR = length(errFCTR{ii});
    errFCTRmat(1:curLenFCTR,ii) = errFCTR{ii};
    curLenCDWS = length(errCDWS{ii});
    errCDWSmat(1:curLenCDWS,ii) = errCDWS{ii};
    curLenCCWS = length(errCCWS{ii});
    errCCWSmat(1:curLenCCWS,ii) = errCCWS{ii};
    curLenSPC = length(errSPC{ii});
    errSPCmat(1:curLenSPC,ii) = errSPC{ii};
end
%% Plot with means
shift=0:0.05:3;
figure, hold on
plot(shift',nanmedian(errFCTRmat),'b') 
plot(shift',nanmedian(errCDWSmat),'g') 
plot(shift',nanmedian(errCCWSmat),'m') 
plot(shift',nanmedian(errSPCmat),'r'), title('FCTR vs. CDWS vs. CCWS vs. SPC')
%% save
save
%% boxplot between CDWS vs CCWS
figure, boxplot(errCDWSmat,'plotstyle','compact','Color','g')
hold on
boxplot(errCCWSmat,'plotstyle','compact','Color','m')
ylim([0 0.21])
%% boxplot between FCTR vs SPC
figure,boxplot(errFCTRmat,'plotstyle','compact','Color','b')
hold on
boxplot(errSPCmat,'plotstyle','compact','Color','r')
ylim([0 0.21])
%% save
save