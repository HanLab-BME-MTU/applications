function [nsl, status]=trackFrame(mov,idlist,dataProperties)
%TRACKFRAME main function for single frame tracker
%
% SYNOPSIS  [nsl, status]=trackFrame(mov,idlist,direction)
%
% INPUT mov : raw microscopy data (2 frames)
%            sl:    list of spots from spotfinder (2 frames)
%
% OUTPUT nsl : tracked frame spots info
%   
% c: 21/1/02	dT

% Constants
%global PIXELSIZE_XY PIXELSIZE_Z FT_SIGMA MAXSPOTS;

if nargin <3 | isempty(dataProperties)
    %load_sp_constants;
    error('Missing dataProperties');
else
    
    PIXELSIZE_XY = dataProperties.PIXELSIZE_XY;
    PIXELSIZE_Z = dataProperties.PIXELSIZE_Z;
    FT_SIGMA = dataProperties.FT_SIGMA;
    MAXSPOTS = dataProperties.MAXSPOTS;
end


PRECISION=[0.01 0.01 0.01];
MAXITER=400;
NUM_PARMS=3;        % currently only translation in 3D
INITIALSCALING = 0.9; % inital scaling for displacement of fused spots


%debug constants
DEBUG=0;
LTYPE={'-','-.',':','--'};
WHATPARMS=[3];

optimmethod=1;  % 1= lsq, 2=lsqnonlin 

status.iterCt =1;
status.msg='Ok.';
% swi=0;

%prefilter frames
%fmov=filtermovie(mov,FILTERPRM);

% init vars
imgTemplate=mov(:,:,:,1,1);
imgTarget=mov(:,:,:,1,2);


%resize data if necessary (odd size)
idx=~rem(size(imgTarget),2);
imgTarget=imgTarget(1:end-idx(1),1:end-idx(2),1:end-idx(3));
imgTemplate=imgTemplate(1:end-idx(1),1:end-idx(2),1:end-idx(3));

%FILTER TARGET TEST
%imgTemplate=filtermovie(imgTemplate,[0.5 0.5 0.5 3 3 3]);
%imgTarget=filtermovie(imgTarget,[0.5 0.5 0.5 3 3 3]);


fSze=size(imgTarget);
centerPos=ceil(fSze/2);


%pixel 2 micron
p2m=ones(MAXSPOTS,1)*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];

%number of tags from spfinder
nspCur=max(idlist(1).linklist(:,2));
nspNext=max(idlist(2).linklist(:,2));

%sort for tag colors
idlist(1).linklist=sortrows(idlist(1).linklist,4);
idlist(2).linklist=sortrows(idlist(2).linklist,4);

%only different spots in source
[tval,uniqSpotIdx,invUniqSpotIdx] = unique(idlist(1).linklist(:,2));

nsl=idlist(2);
%[nsl.spot(1:nspNext).mult]=deal(-1);


%coordinate in mu (conver to matlab convention x<->y)
coordsCurMu=idlist(1).linklist(uniqSpotIdx,9:11);
coordsCurMu=[coordsCurMu(:,2) coordsCurMu(:,1) coordsCurMu(:,3)];

%use corresponding coords
coordsNextMu=idlist(2).linklist(uniqSpotIdx,9:11);
coordsNextMu=[coordsNextMu(:,2) coordsNextMu(:,1) coordsNextMu(:,3)];

%create a scaling vector for the initial displacement parameters (=1 for isolated tags, =INITIALSCALING for fused tags)
initialParamScaling=ones(1,length(uniqSpotIdx));
[sps sIdx]=sort(idlist(2).linklist(uniqSpotIdx,2));
dsps=diff(sps);
idsps=find(dsps==0);
if ~isempty(idsps)
    initialParamScaling(sIdx(idsps))=INITIALSCALING;
    initialParamScaling(sIdx(idsps+1))=INITIALSCALING;
end
%convert coordinates to pixels 
coordsCurPix=coordsCurMu./p2m(1:nspCur,:);
coordsNextPix=coordsNextMu./p2m(1:nspCur,:);

dnsp=nspCur - nspNext;

%create masks and perform spot segmentation
match=maskSpots(idlist,fSze,dataProperties);

if DEBUG
    global path;
    path=[];
end;

% prepare match spots & init params
for s=1:nspCur 
    match(s).fSze=fSze;
    % currently 2D noise estimation from central section
    match(s).sigmaNoise=imNoiseEstim(imgTemplate(:,:,centerPos(3),:,1));
    if DEBUG
        ['a priori noise estimation: ' num2str(match(s).sigmaNoise)];
    end;
    match(s).center=coordsCurPix(s,:);
    %match(s).amp=idlist(1).linklist(s,8);
    %    match(s).data=(imgTemplate(sub2ind(fSze,match(s).coords(:,1),match(s).coords(:,2),match(s).coords(:,3)))-sl(1).sp(s).bg).*match(s).ratio;
    match(s).data=(imgTemplate(sub2ind(fSze,match(s).coords(:,1),match(s).coords(:,2),match(s).coords(:,3)))-mean(imgTemplate(:))).*match(s).ratio;
    %    match(s).amp=max(match(s).data);
    % currently just translation
    match(s).parms=initialParamScaling(s)*[coordsNextPix(s,:)-coordsCurPix(s,:)];
    % Set intial parameter change
    match(s).dpar=2*PRECISION;
    
    osc(s).mask=ones(1,NUM_PARMS);

    
    %init osc_check
    osc(s).dp=ones(3,NUM_PARMS);
    
    % FILTER MATCH
    %dermine minimal surrounding box
    %      minV=[min(match(s).coords(:,1)) min(match(s).coords(:,2)) min(match(s).coords(:,3))]-1;
    %      maxV=[max(match(s).coords(:,1)) max(match(s).coords(:,2)) max(match(s).coords(:,3))]+1;
    %      gradSize=maxV-minV+1;
    %      
    %     
    %     %shift coord
    %      shiftC=([minV(1)-1, minV(2)-1, minV(3)-1]);
    %      gcords=match(s).coords-ones(size(match(s).coords,1),1)*shiftC;
    %      mImg=zeros(gradSize);
    %      mImg(sub2ind(gradSize,gcords(:,1),gcords(:,2),gcords(:,3)))=match(s).data;
    %      fmImg=fastGauss3d(mImg,[1 1 1],[5 5 5]);
    %      match(s).data=fmImg(sub2ind(gradSize,gcords(:,1),gcords(:,2),gcords(:,3)));
%    lastXi(s)=10;
end;


%--------------------main tracking loop start-----------------------
iterCt=0;

if optimmethod==1
    if DEBUG
        noH=figure('Name','chi squared');
        hold on;
        dpH=figure('Name','parm change');
        hold on;
        parmH=figure('Name','parm value');
        hold on;
    end;
    %----------------std LSQ
    
    %Counter for iteration
    %h=1;
    PREMAT=PRECISION'*ones(1,length(match));
    PREMAT=PREMAT(:)';
    while( (any(abs([match.dpar].*[osc.mask])>PREMAT) & iterCt<MAXITER & status.iterCt==1))
        %update counter redraw figures
        iterCt=iterCt+1;
        drawnow;
        %Compute new coords according to model
        model=transSpot(match,dataProperties);
        for s=1:length(model)
            %check coordinates
            if any((model(s).coords(:,1)<1) | (model(s).coords(:,2)<1) | (model(s).coords(:,3)<1) | (model(s).coords(:,1)>fSze(1)) | (model(s).coords(:,2)>fSze(2))| (model(s).coords(:,3)>fSze(3)))
                status.msg='Coordinates out of target image';
                %status=2;
                model(s).coords(:,1)=max(model(s).coords(:,1),1);
                model(s).coords(:,1)=min(model(s).coords(:,1),fSze(1));
                model(s).coords(:,2)=max(model(s).coords(:,2),1);
                model(s).coords(:,2)=min(model(s).coords(:,2),fSze(2));
                model(s).coords(:,3)=max(model(s).coords(:,3),1);
                model(s).coords(:,3)=min(model(s).coords(:,3),fSze(3));               
            end
            % SOMETHING LIKE THIS LINE FOR ESTIMATED BACKGROUND: model(s).data=interpn(imgTarget-sl(2).sp(idmap(1).sp{s}).bg,model(s).coords(:,1),model(s).coords(:,2),model(s).coords(:,3)).*model(s).ratio;
            model(s).data=interpn(imgTarget-mean(imgTarget(:)),model(s).coords(:,1),model(s).coords(:,2),model(s).coords(:,3),'cubic').*model(s).ratio;
            
            %compute gradient
            % FILTER MODEL
            %dermine minimal surrounding box
            minV=[min(match(s).coords(:,1)) min(match(s).coords(:,2)) min(match(s).coords(:,3))]-1;
            maxV=[max(match(s).coords(:,1)) max(match(s).coords(:,2)) max(match(s).coords(:,3))]+1;
            gradSize=maxV-minV+1;
            
            
            %shift coord for gradient calculation
            shiftC=([minV(1)-1, minV(2)-1, minV(3)-1]);
            gcords=match(s).coords-ones(size(match(s).coords,1),1)*shiftC;
            
            tImg=zeros(gradSize);
            tImg(sub2ind(gradSize,gcords(:,1),gcords(:,2),gcords(:,3)))=model(s).data;
            %FILTER
            %ftImg=fastGauss3d(tImg,[1 1 1 ],[5 5 5]);            
            %model(s).data=ftImg(sub2ind(gradSize,gcords(:,1),gcords(:,2),gcords(:,3)));
            %END FILTER
            model(s).residual=(match(s).data -model(s).data);
            
            %tImg=filtermovie(tImg,[1 1 0.7 3 3 3]);
            [FX,FY,FZ] = gradient(tImg);
            FX=FX(sub2ind(gradSize,gcords(:,1),gcords(:,2),gcords(:,3)));
            FY=FY(sub2ind(gradSize,gcords(:,1),gcords(:,2),gcords(:,3)));
            FZ=FZ(sub2ind(gradSize,gcords(:,1),gcords(:,2),gcords(:,3)));
            
            model(s).iGrad=[FY(:) FX(:) FZ(:)];
            
            %current estimated aposteriori noise
            model(s).sigma0=sqrt(sum((model(s).residual).^2./(model(s).ratio.^2+match(s).ratio.^2))/(length(model(s).residual)-NUM_PARMS));
            if DEBUG
                %['a posteriori noise estimation of tag ' num2str(s) ':  ' num2str(model(s).sigma0)];
            end;
        end;
        % error with current params: e = mA*params-mB
        %dp=blkdiag(model(1).iGrad, model(2).iGrad)\blkdiag(model(1).residual, model(2).residual);
        %dp
        
        %        msk(1,:)=[abs(match(1).dpar(1))>PRECISION(1) abs(match(s).dpar(2))>PRECISION(2) abs(match(s).dpar(3))>PRECISION(3)];
        
        %         if swi==1 | all(abs(match(1).dpar)<PRECISION)
        %             swi=1;
        %             msk(1,:)=[0 0 0];
        %             msk(2,:)=[1 1 1];
        %         end;
        %         if all(abs(match(2).dpar)<PRECISION)
        %             swi=0;
        %             msk(1,:)=[1 1 1];
        %             msk(2,:)=[0 0 0];
        %         end;
        
        for s=1:length(model)
            match(s).dpar= model(s).iGrad\model(s).residual;
                        
            match(s).dpar=match(s).dpar';
                        
%            curXi(s)=sum(model(s).residual.^2)/length(model(s).residual);
%            msign(s)=(curXi(s)<lastXi(s))-(curXi(s)>=lastXi(s));
            match(s).parms=match(s).parms+(match(s).dpar.*osc(s).mask);
%            lastXi(s)=curXi(s);
            
            %check for oscillation 
            %osc(s).dp=[osc(s).dp(2:2,:); match(s).dpar];
            %osc(s).mask=abs(osc(s).dp(2:2,:))<=abs(osc(s).dp(1:1,:));
            %osc(s).mask=prod(osc(s).mask,1);
            
            if DEBUG
                %for DEBUG
                debInfo(s,iterCt).stno=(sum(model(s).residual.^2)/length(model(s).residual));
                
                match(s).dpar;
                match(s).center+match(s).parms;
                norm(model(s).residual);
                debInfo(s,iterCt).std=std(model(s).residual);
                debInfo(s,iterCt).stno;
                debInfo(s,iterCt).dp=match(s).dpar;
                debInfo(s,iterCt).mpar=match(s).parms+match(s).center;
                dp= cat(1,debInfo(s,:).dp);
                mpar=cat(1,debInfo(s,:).mpar);
                
                for i=WHATPARMS
                    figure(dpH);
                    plot(dp(:,i),LTYPE{s},'Color',[1/i .2 .2]);
                    figure(parmH);
                    plot(mpar(:,i),LTYPE{s},'Color',[1/i .2 .2]);
                end;
            end;
        end;

        if DEBUG
            figure(noH);
            plot([debInfo(1,:).stno],'r');
            hold on;
            plot([debInfo(2,:).stno],'-.');
            plot(([debInfo(1,:).stno]+[debInfo(2,:).stno])/2,'g');
            drawnow;
        end;
        
    end;
    
else
    % 
    %-------------- LSQNONLIN --------------------------------
    % init lsqnonlin
    iterCt=0;
    options = optimset('Jacobian','off','Display','on','LevenbergMarquardt','on','LargeScale','on');
    lb=[-9 -9 -9];
    ub=[9 9 9];
    %loop until parms change < PREC
    %   while( (any(abs([match.dpar])>PRECISION) & iterCt<MAXITER) )
    iterCt=iterCt+1;
    %            [parms,resnorm] = lsqnonlin(@modelError,[match.parms],lb,ub,options,match,s,imgTarget,sl,idmap);
    [parms,resnorm] = lsqnonlin(@modelError,[match.parms],[],[],options,match,imgTarget,dataProperties);
    for s=1:length(match)          
        match(s).dpar=parms((s-1)*3+1:s*3)-match(s).parms;
        match(s).dpar; %for DEBUG
        match(s).parms=parms((s-1)*3+1:s*3);
        
        debInfo(s,iterCt).stno=resnorm;
    end;
    model=transSpot(match,dataProperties);
    for s=1:length(match)          
        model(s).data=interpn(imgTarget-mean(imgTarget(:)),model(s).coords(:,1),model(s).coords(:,2),model(s).coords(:,3),'cubic').*model(s).ratio;
        model(s).residual=(match(s).data -model(s).data);
        model(s).sigma0=sqrt(sum((model(s).residual).^2./(model(s).ratio.^2+match(s).ratio.^2))/(length(model(s).residual)-NUM_PARMS));
    end
    
    if DEBUG
        %iterCt
        plot([debInfo(1,:).stno],'r');
        hold on;
        %plot([debInfo(2,:).stno],'-.');
    end;
end;
%--------------------main tracking loop end-----------------------

% for s=1:length(match);
%     % reject if significant (apriori/aposteriori) noise differences
%     if  ~acceptMatch(model(s).sigma0,match(s).sigmaNoise,length(model(s).residual)-NUM_PARMS)
%         status.iterCt=-1;
%         status.msg='tracking rejected';
%     end;
% end;

%check for max iteration
if iterCt<MAXITER
    %check for other error
    if status.iterCt>0
        if optimmethod==1
            %refill nsl
            for s=1:length(match);
                model(s).QMatrix=(model(s).iGrad'*model(s).iGrad)^-1;
            end;
        end;
        %check for significance
        %        matchOK=acceptTrackedSpotDist(match,model);
        for s = 1:length(invUniqSpotIdx);
            %             if matchOK(invUniqSpotIdx(s))
            tcoord= match(invUniqSpotIdx(s)).center+match(invUniqSpotIdx(s)).parms;
            %             else
            %                 tcoord=coordsNextPix(invUniqSpotIdx(s),:);
            %             end;
            tcoord=[tcoord(2) tcoord(1) tcoord(3)].*p2m(1,:);
            nsl.linklist(s,9:11)=tcoord;
            if optimmethod==1
                nsl.info.trackQ_Pix=blkdiag(model.QMatrix);
            end;
            if DEBUG
                nsl.linklist(s,:);   % for DEBUG
            end;
        end;
        
        spmax=1;
        for s=1:size(nsl.linklist,1)
            %test for spot separation
            fusionIdx=find(all(ones(size(nsl.linklist,1),1)*nsl.linklist(s,9:11)==nsl.linklist(:,9:11),2));
            %update spot color
            nsl.linklist(s,3)=sum(nsl.linklist(fusionIdx,4));
            if s==fusionIdx(1)
                % and spot number of second spot
                nsl.linklist(s,2)=spmax;
                spmax=spmax+1;
            else
                nsl.linklist(s,2)=nsl.linklist(fusionIdx(1),2);
            end;
        end;
    end;
else
    status.msg='Max iterations exceeded.';
end;
status.iterCt=iterCt;