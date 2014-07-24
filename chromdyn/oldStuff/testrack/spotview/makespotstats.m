function makespotstats
% MAKESPOTSTATS evalutes and displays slist content

%CONST DEFINITIONS
load_sp_constants

quit=0;
hold on;
LINESTYLE={'-','--',':'};
lc=1;

while (~quit)
    nsp=[];
    meansnr=[];
    minsnr=[];    
    lmdist=[];
    mmdist=[];
    lm=0;
    mm=0;
    
    spstat=zeros(1,10);
    
    % read data file
    [fname,path]=uigetfile('','select result file');
    if (isempty(fname))
        return;
    end;
    cd(path);
    filename=[path fname];
    load(filename);
    
    nframes=length(slist);
    for t=1:nframes
        spb=[];
        
        %count spots
        if (~isempty(slist(t).sp))
            nsp(t)=size(cat(1,slist(t).sp.cord),1);
        else
            nsp(t)=0;
        end;
        spstat(nsp(t)+1)=spstat(nsp(t)+1)+1;
        
        %loopt through spots
        for s=1:nsp(t)
            %assign spb coords
            if strcmp(slist(t).sp(s).type,'spb');
                spb=[spb; slist(t).sp(s).cord];
            end;
            %separate locmac (lm) spots and multi-mix (mm) spots
            if slist(t).sp(s).mult
                mm=mm+1;
            else
                lm=lm+1;
            end;
            if any([slist(t).sp.mult])
                mmdist=[mmdist min(nonzeros(slist(t).distMat))];
            else
                lmdist=[lmdist min(nonzeros(slist(t).distMat))];
            end;
        end;
        
        %fill snr
        if(~isempty(slist(t).statistics))
            meansnr(t)=mean(slist(t).statistics.snr);
            minsnr(t)=min(slist(t).statistics.snr);
        else
            if nsp(t)==4
                meansnr(t)=meansnr(t-1);
                minsnr(t)=minsnr(t-1);
            else
                meansnr(t)=NaN;
                minsnr(t)=NaN;
            end
        end;
        
        %comp pol2poldist
        %         if size(spb,1)==2
        %             dif=(spb(1,:)-spb(2,:)).*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];
        %             pol2poldist(t)=sqrt(sum(dif.^2));
        %         else
        %             if t>1
        %                 pol2poldist(t)=pol2poldist(t-1);
        %             else
        %                 pol2poldist(t)=0;
        %             end;
        %         end;
        
        if( ~isempty(slist(t).distMat) & ~all(slist(t).distMat==0))
            pol2poldist(t)=max(max(slist(t).distMat));
        else
            pol2poldist(t)=pol2poldist(t-1);
        end;
        
        
    end;
    %norm
    tot=sum(nsp);
    mm=mm/tot;
    lm=lm/tot;
    spstat=spstat/nframes;
    
    
    %output
    disp(['number of time points in ',path,': ',num2str(nframes)]);
    
    disp(sprintf('%s  %s','# spots','%of total time points'));
    disp(sprintf('%i            %1.2f\n',[0:9'; spstat]));
    
    disp(sprintf('%s%1.5f','%of mix model spots: ',mm));
    if ~isempty(lmdist)
        disp(sprintf('%s%1.5f, %1.5f','meanmin and minmin dist of loc max spots: ',mean(lmdist),min(lmdist)));
    end;
    if ~isempty(mmdist)
        disp(sprintf('%s%1.5f, %1.5f','meanmin minmin dist of mix model spots: ',mean(mmdist),min(mmdist)));
    end;
    
    %plot
    hp1=line(1:nframes,nsp,'Color','b');
    ax1 = gca;
    axis([0 100 0 max(nsp)+1]);
    set(ax1,'FontSize',20,'YColor','b');
    box off;
    %set(hp1,'LineStyle',LINESTYLE{mod(lc-1,length(LINESTYLE))+1},'Color','b');
    ax2 = axes('Position',get(ax1,'Position'),...
        'FontSize',20,...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor','r');
    hp2=line(1:nframes,meansnr,'Color','r','LineStyle','--');
    hp3=line(1:nframes,minsnr,'Color','r','LineStyle',':');
    %set(hp2,'LineStyle',LINESTYLE{mod(lc-1,length(LINESTYLE))+2},'Color','g');
    axis([0 100 0 max(meansnr)]);
    box off;
    ppfH=figure;
    plot(pol2poldist)
    hold on;
    brob=robustfit(1:nframes,pol2poldist);
    plot(1:nframes,brob(1)+brob(2)*[1:nframes],'r-')
    lc=lc+1;
    
    %     button = questdlg('Do you want to add a data set?',...
    %         'Continue Operation','Yes','No','No');
    %     if strcmp(button,'No')
    quit=1;
    %     end
end;