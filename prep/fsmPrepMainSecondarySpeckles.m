function [IMfinal,candsTot]=fsmPrepMainSecondarySpeckles(I,strg,counter,noiseParam,Speckles,enhTriang,fsmParam,oI)

% fsmPrepMainSecondarySpeckles is the main function of the fsmPrepSecondarySpeckles sub-module
%
% Plots PSFs on the positions of the primary speckles (found by locmax operator) and substructs
% them from the filtered data. Appled again on the resulting image, the locmax-operator finds 
% new (secondary) speckles (intensity and distance significance tests applied)
%
%
% SYNOPSIS   [IMfinal,candsTot]=fmsPrepMainSecondarySpeckles(I,strg,counter,noiseParam,Speckles,enhTriang)
%
% INPUT      I          :  filtered image
%            strg       :  format string for the correct file numbering
%            counter    :  image number
%            noiseParam :  noise parameters for statistical speckle selection
%            Speckles   :  (1) contains information for the hierarchical level 
%                          (2) minimal increase in (%) of new speckles
%                          (before stopping)
%            enhTriang  :  turns on enhanced triangulation for Matlab Version < 6.5
%            fsmParam   :  (optional) fsmParam structure for SpeckTackle
%            oI         :  original image necessary for subpixel accuracy
%                           determination; is normalized, but NOT filtered
%
% OUTPUT     candsTot   :  augmented cands structure (see fsmPrepTestLocalMaxima.m)
%            IMfinal    :  local maxima map
%
%
%
% DEPENDENCES   fsmPrepMainSecondarySpeckles uses { fsmPrepConfirmSpeckles, fsmPrepSubstructMaxima, 
%                                         fsmPrepCheckDistance, fsmPrepUpdateImax, fsmPrepCheckInfo}
%               fsmPrepMainSecondarySpeckles is used by { fsmPrepMain }
%
% Alexandre Matov, November 7th, 2002

if nargin==0
    locDEBUG=1;
    DEBUG=0;
    Speckles=[3 0.0/100];
    %     Speckles=[8 0.0/100]; % default
    % load first image
    [fileName,dirName] = uigetfile('*.tif','Choose an image');
    I=imread([dirName,filesep,fileName]);
    SIG=1.88;
    enhTriang=0;
    %     I=double(I); %??
    IG=prepareRowData(I,SIG);
    strg=[];
    counter=0;
    shift=4;%15
    
    % format: [Q/GaussRatio sDN beta I0 Q]
    %      noiseParam=[1.96/1.12985 0.00028 1e-4 0.00495 1.96];  % WT 2s 14bit GOOD ONE (but run with division like 8bit!)
%          noiseParam=[1.96/1.17700 0.00039 2e-4 0.00531 1.96];  % WT 2s 14bit EVEN BETTER (but run with division like 8bit!)
%          noiseParam=[1.96/1.20461 0.01695 1e-5 0.30953 1.96];  % meta spindle WT 2s 8bit ????
    %      noiseParam=[1.96/1.21230 0.01036 1e-4 0.27200 1.96];  % meta spindle WT 2s 8bit BLACK BG
    
%          noiseParam=[1.96/1.25689 0.01200 1e-4 0.29697 1.96];  % meta spindle WT 4s 8bit GOOD ONE
    
%          noiseParam=[1.96/3.55 0.03656 1e-4 0.14540 1.96]; % actin retrograde flow
    
%         noiseParam=[1.96/3.55 5.4103e-004 2e-4 0.0295 1.96];% actin 158a 14bit
    % noiseParam=[1.96/3.55 5.4103e-004 2e-4 0.0295 1.96];% actin 07 14bit NOISE PARAMS NOT GOOD
    
    %     noiseParam=[1.96/2 0.0002 1e-4 0.003306872 1.96]; % image Michael AOTF + EPI (naglaseni)
    %     noiseParam=[1.96/1.03897826 0.00059316 1e-4 0.03213289 1.96]; % getmodulation UNIVERSAL
    noiseParam=[1.4700    0.0001    0.0002    0.0313    1.9600]; % 488 S1 AOTF488_15_1.tif
%     noiseParam=[1.96/1.0993 0.0002 1.e-4  0.0316 1.96]; % 488 S1 AOTF488_20_1.tif
%     noiseParam=[1.96/1.0597 0.0004 1.e-4  0.0318 1.96]; % 488 S1 AOTF488_25_1.tif
%     noiseParam=[1.96/1.0485 0.0004 1.e-4  0.0320 1.96]; % 488 S1 AOTF488_32_1.tif
%       noiseParam=[1.96/1.7734 0.0008 7.e-4  0.0311 1.96]; % 488 S1 Epi488_1.tif



% noiseParam=[1.1660    0.0006    0.0009    0.0304    1.9600]; % 488 STACK
      
    %     noiseParam=[1.96/1.03897826 0.00059316 1e-4 0.03213289 1.96]; % image Michael AOTF
    
    %     noiseParam=[1.96/2.02281245 0.00074871 7e-4 0.03086343 1.96]; % image Michael EPI 1
    % noiseParam=[1.96/2.26574167 0.00070565 13e-4 0.03062333 1.96]; % image Michael EPI 2
    
%         noiseParam=[1.96/1 0.02 2e-4 0.0981 1.96]; % WT 10s 8bit RED
    % change sigma dark noise
    % noiseParam(2)=noiseParam(2)/5.2540;
    
    %     noiseParam=[1.96/2.1090 0.02085 1e-4 0.09317 1.96];% Confocal images1 almost o.k. 
    %     noiseParam=[1.96/2.75734 0.01971 1e-4 0.08296 1.96];% Confocal images2 too weak
    %     noiseParam=[1.96/1.22425 0.04055 1e-4 0.15979 1.96];% Confocal images3
    SAVEINFO=0;
else
    IG=I;
    SAVEINFO=1;
    if strg == 0 % if you provide all the fields but dont wanna write to disc
        SAVEINFO=0;
    end
    locDEBUG=0;
    DEBUG=0;
end

if nargin==6
    fsmParam=[];
    userROIbw=[];
end

if ~isempty(fsmParam)
    if fsmParam.prep.drawROI~=0 % Either 1 (drawn) or 2 (loaded)
        
        % Load user-defined ROI from disk
        ROIname=[fsmParam.main.path,filesep,'userROI.mat'];
        if exist(ROIname)==2 % File found
            tmp=load(ROIname);
            userROIbw=tmp.userROIbw;
            clear tmp;
        else
            userROIbw=[];
        end
    else
        userROIbw=[];    
    end
end

% SIG=1.60;
% SIG=1.77;
SIG=1.88; % for the twice convolved image (or 1.77)

% local minima
Imin=locmin2d(IG,[3,3]);

% intial (filtered) image
[yi,xi,y,x,Imax,candsP,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam,enhTriang,userROIbw); % TO DO: update cands

aux=length(candsP);
for i=1:aux
    candsP(i).speckleType=1;
end

candsTot=candsP;

Inew=IG;
candsS=candsP;
HierLevel=2;

while HierLevel<=Speckles(1) & length(candsS)>(Speckles(2)*length(candsTot)) & length(find([candsS.status]==1))>0
    
    [Inew,Imaxima,nmB]=fsmPrepSubstructMaxima(Inew,Imax,SIG,candsS); % prednite Cands
    [yni,xni,yn,xn,Imax,candsS]=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,enhTriang,triMin,pMin,IG,userROIbw);
    %     [yni,xni,yn,xn,Imax,candsS]=fsmPrepConfirmSpeckles(Inew,Imin,noiseParam,enhTriang);
    
    aux=length(candsS);
    for i=1:aux
        candsS(i).speckleType=HierLevel; % type flag
    end
    
    candsS=fsmPrepCheckDistance(candsS,candsTot); 
    
    HierLevel=HierLevel+1;
    
    if ~isempty(candsS)
        candsTot=cat(2,candsTot,candsS); % concatenating secondary and primary cands structures
    end
    
end

% remove repetitions because of secondary speckles apearing on the same positions as primary (because of floating background)
candsTot=fsmPrepCheckInfo(candsTot);

% obtain updated IM from candstTot
IMfinal=zeros(size(IG));
for i=1:length(candsTot)
    if candsTot(i).status==1
        IMfinal(candsTot(i).Lmax(1),candsTot(i).Lmax(2))=candsTot(i).ILmax;
    end
end
[yMfinal,xMfinal]=find(ne(IMfinal,0));


% Save speckle information (cands and locMax) to disk for future use
if SAVEINFO==1
    locMax=IMfinal;
    cands=candsTot;   
    indxStr=sprintf(strg,counter);
    eval(strcat('save cands',filesep,'cands',indxStr,'.mat cands;')); % Save speckle info
    eval(strcat('save locMax',filesep,'locMax',indxStr,'.mat locMax;')); % Save loc max positions
    
end

%-----------------------------------------
% Estimate subpixel positions of speckles
%-----------------------------------------
% The size of the GaussKernel used for filtering the image (in
% fsmPrepareImage) has, as of March 2005, been set to the real value of 
% psfsigma as determined by the physical parameters,
% psfsigma=0.21*(lambda/NA)/pixelsize
% Thus, the the mixture model fitting may be performed without loss of
% information on the filtered image (called I), rather than the original
% image (called oI)
%
% However, it is important to note that performing the Gauss fit in the
% mixture model on the filtered image (which is broadened), rather than on
% the original one, requires to modify the sigma of the mixture-model fit!


if fsmParam.prep.subpixel==1
    cands=candsTot; 
    psfsigma     = fsmParam.prep.psfSigma;       % true physical sigma of the image point-spread function, caluclated by sigma=0.21*(lambda/NA)/pixelsize
    filtersigma  = fsmParam.prep.filterSigma;    % sigma used for the low-pass filtering; except where specifically
                                                % stated differently by the user, filtersigma should have the same value as psfsigma; 
                                                % for filtersigma>psfsigma, image information is lost during filtering!!                                            % same value as 
    %mixture model Gauss sigma (mmsigma) is calculated from psfsigma and
    %filtersigma; in the usual case where psfsigma=filtersigma, then
    %mmsigma=sqrt(2)*psfsigma
    mmsigma=sqrt(psfsigma^2+filtersigma^2);
    image=I;
    disp(['psfsigma=',num2str(psfsigma),'   filtersigma=',num2str(filtersigma),'   mixmodsigma=',num2str(mmsigma)]);
    disp('calculating sub-pixel locations...');
    [candsSP] = candsToSubpixelN(image,cands,mmsigma);
    eval( (strcat('save cands',filesep,'cands',indxStr,'_spa.mat candsSP;')) );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG FIGURES:

if locDEBUG==1
    if DEBUG==1
        
        %         length(find(IM~=0))  
        length(candsTot) 
        
        c1=0;
        for i=1:length(candsTot)
            if candsTot(i).status==1
                c1=c1+1;
                SNR(c1)=candsTot(i).deltaI/candsTot(i).sigmaLmax;
            end
        end
        figure(100),plot(SNR);
        MEAN_SNR=mean(SNR)
        
        c2=0;
        for j=1:length(candsTot)
            if candsTot(j).status==0
                c2=c2+1;
                SNRins(c2)=candsTot(j).deltaI/candsTot(j).sigmaLmax;
            end
        end
        figure(101),plot(SNRins);
        MEAN_SNR_ins=mean(SNRins)
        
        PRIMARY=size(y,1)
        Primary_rejected_noiseModel=size(yi,1)-size(y,1)
        PRIMARY_rejected_procentage=Primary_rejected_noiseModel/size(yi,1)*100
        
        % the figure is showing all local maxima and also only the significant ones
        figure(1);
        imshow(IG(1+shift:end-shift,1+shift:end-shift),[]);
        hold on;
        plot(xi-shift,yi-shift,'y*');
        plot(x-shift,y-shift,'r.');
        hold off;
        title('all the local maxima in yellow and the significant ones in red asteriques');
        
        % local maxima points in the raw data
        figure(2);
        surf(Imax);
        title('poinst at positions of the local maxima but with intensity DELTA.I');
        
        figure(19);
        imshow(IG(1+shift:end-shift,1+shift:end-shift),[]);
        hold on;
        plot(x-shift,y-shift,'r*');
        plot(xn-shift,yn-shift,'g*');
        hold off;
        title('red asteriques are the initial local maxima and the green dots are the secondary speckles found');
        
        figure(20);
        imshow(IG(1+shift:end-shift,1+shift:end-shift),[]);
        hold on;
        plot(x-shift,y-shift,'r*');
        hold off;
        title('red asteriques are the initial local maxima and the green and blue dots are the secondary speckles found');
        
        figure(21);
        imshow(IG(1+shift:end-shift,1+shift:end-shift),[]);
        hold on;
        plot(x-shift,y-shift,'r*');
        hold off;
        title('red asteriques are the initial local maxima and the green dots are the secondary speckles found');
        
        figure(22);
        imshow(IG(1+shift:end-shift,1+shift:end-shift),[]);
        hold on;
        
        for i=1:length(candsTot)
            if candsTot(i).status==1
                plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'y.');
            end
        end
        hold off;
        title('red asteriques - IM; yellow dots - cands');
    end
    
%     figure
%     imshow(IG(1+shift:end-shift,1+shift:end-shift),[]);
%     hold on
%     for i=1:length(candsTot)
%         if candsTot(i).status==1
%             plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'y.');
%         end
%     end
%     hold off
%     title('all speckles yellow');
    % END WAS HERE
    %end
    
    SECONDARY=0;
    PRIMARY=0;
    TERTIARY=0;
    QUATTRO=0;
    CINQUE=0;
    SHEST=0;
    SEDEM=0;
    OSEM=0;
    shift=4;%15
    
    
    figure,imshow(IG(1+shift:end-shift,1+shift:end-shift),[]);
    hold on
    for i=1:length(candsTot)
        if candsTot(i).status==1
            switch candsTot(i).speckleType
                case 1 
                    plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'g.'); PRIMARY=PRIMARY+1;
                case 2 
                    plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'r.'); SECONDARY=SECONDARY+1;
                case 3 
                    plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'b.'); TERTIARY=TERTIARY+1;
                case 4 
                    plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'c*'); QUATTRO=QUATTRO+1;
                case 5
                    plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'m*'); CINQUE=CINQUE+1;
                case 6
                    plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'y.'); SHEST=SHEST+1;
                case 7
                    plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'y*'); SEDEM=SEDEM+1;
                case 8
                    plot(candsTot(i).Lmax(2)-shift,candsTot(i).Lmax(1)-shift,'y*'); OSEM=OSEM+1;
                otherwise
                    error('wrong classification of secondary speckles');
            end
        end
    end
    hold off
    title('the green spots denote significant speckles');
        title('green-primary red-secondary blue-tertiary cyan-quattro magenta-cinque yellow-shest');
    
        PRIMARY
        SECONDARY
        TERTIARY
        QUATTRO
        CINQUE
        SHEST
        SEDEM
        OSEM
    
    % MOVE IT BACK DOWN
end