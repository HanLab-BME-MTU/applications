function monitorStat(stat, opt)

load_sp_constants;
global FT_SIGMA PIXELSIZE_XY PIXELSIZE_Z;

if any(opt==1)
    figure;
    hold on;
    curcol=1;
    for nsct = 1:size(stat,1);
        ct=1;
        deterr=[];
        trackerr=[];
        odist=[];
        for i=1:size(stat,2)
            if stat(nsct,i).foundStat>0.2
                %orig=[stat(nsct,i).origSp2(2) stat(nsct,i).origSp2(1)  stat(nsct,i).origSp2(3)];
                %perr=[cat(1,stat(nsct,i).sp2.cord)-ones(length(stat(nsct,i).sp2),1)*orig].*(ones(length(stat(nsct,i).sp2),1)*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z]);
                %sqperr=sum(perr.^2,2);
                deterr(ct)=1000*sqrt(sum(stat(nsct,i).detecter1)/length(nonzeros(stat(nsct,i).detecter1)));
                trackerr(ct)=1000*sqrt(sum(stat(nsct,i).tracker1)/length(nonzeros(stat(nsct,i).tracker1)));
                odist(ct)=stat(nsct,i).origDist/(FT_SIGMA(2)*PIXELSIZE_XY);
                ct=ct+1;
            end;
        end;
        
        plot(odist,deterr,'Color',[curcol .2 .2]);
        plot(odist,trackerr,'-.','Color',[.2 .2 curcol]);
        curcol=curcol-1/size(stat,1);
        pause;
    end;
    box;
    
end;



if any(opt==2)
    
    dist=1;
    figure;
    hold on;
    deterr=[];
    trackerr=[];
    odist=[];
    for nsct = 1:size(stat,1);
        
        deterr(nsct)=1000*sqrt(sum(stat(nsct,dist).detecter1)/length(nonzeros(stat(nsct,dist).detecter1)));
        trackerr(nsct)=1000*sqrt(sum(stat(nsct,dist).tracker1)/length(nonzeros(stat(nsct,dist).tracker1)));
        snrl(nsct)=stat(nsct,dist).snr;
        
    end;
    plot(snrl,deterr,'r');
    plot(snrl,trackerr,'b-.');
    box;
end;


return;

%------------------------------------------------------

%PLOT distance
figure;
plot([1 length(fdStatistics(noiseCt,distCt).trackfoundDist)],[fdStatistics(noiseCt,distCt).origDist fdStatistics(noiseCt,distCt).origDist]);
hold on;
plot([1:length(fdStatistics(noiseCt,distCt).detectfoundDist)],fdStatistics(noiseCt,distCt).detectfoundDist,'r');
plot([1:length(fdStatistics(noiseCt,distCt).trackfoundDist)],fdStatistics(noiseCt,distCt).trackfoundDist,'b');

%PLOT errors
figure(errDist);
set(errDist,'Name',['SNR=' num2str(snr)]);
deter=[];
tracker=[];
for i = 1:size(fdStatistics,2)
    deter(i)=mean(fdStatistics(noiseCt,i).detecter1);
    tracker(i)=mean(fdStatistics(noiseCt,i).tracker1);
end;
plot([1:size(fdStatistics,2)],deter,'r');
plot([1:size(fdStatistics,2)],tracker,'b');

figure(errSNR);
set(errSNR,'Name',['dist=' num2str(sp_dist)]);
deter=[];
tracker=[];
for i = 1:size(fdStatistics,1)
    deter(i)=mean(fdStatistics(i,distCt).detecter1);
    tracker(i)=mean(fdStatistics(i,distCt).tracker1);
end;
plot([1:size(fdStatistics,1)],deter,'r');
plot([1:size(fdStatistics,1)],tracker,'b');
