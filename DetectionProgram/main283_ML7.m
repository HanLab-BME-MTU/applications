function [] = main283_ML7(iclean,icut,itest);
%reads a whole bunch of original images and for each image it
% Reconstructs a filtered image using iterative filtering from significant
% coefficients ( see Book by Starck et al, p. 66)

%Then calculates a binary image and the product of the binary image and the
%reconstructed image

% It then finds and plots ALL the local maxima on each cluster (at least
% one maximum per cluster)

%If a cluster has more than one maximum, it divides them into primary and
%secondary (The secondaries are listed at the end.)

%It calculates the intensity of each maximum and the total intensity and size (i.e. area) of
%each cluster

%and then creates a movie.

%%========================================================================
%
%   NOTE: this is Dinah's modified version, which runs under
%   more recent versions of matlab and therefore on the LINUX server
%   
%   Last Date of modification: March 14, 2007
%
%=========================================================================

%   INPUT:  iclean (optional) =     input iclean value; if no value is set 
%                                   here, or input is empty, the default 
%                                   value will be 
%                                   iclean=0
%           icut   (optional) =     input icut value; if no value is set 
%                                   here, or input is empty, the default 
%                                   value will be 
%                                   icut=5
%           itest   (optional) =    input icut value; if no value is set 
%                                   here, or input is empty, the default 
%                                   value will be 
%                                   itest=0
% EXAMPLE: main283_ML7(1,[],0); 
%
% NOTE: the default values can be changed below in line 58-60, if desired,
% but generally, different parameter values can be used by entering them in
% the commando line instead of changing the code


close all;
clc;


%=========================================================================
%
%       set default values for parameters
%
%=========================================================================

icleandefault = 0;
icutdefault   = 5;
itestdefault  = 0;


% if iclean is not defined or if its input is [], set iclean to default
if nargin < 1
    iclean = icleandefault;
elseif isempty(iclean)
    iclean = icleandefault;
end


% if icut is not defined or if its input is [], set icut to default
if nargin < 2
    icut = icutdefault;
elseif isempty(icut)
    icut = icutdefault;
end


% if itest is not defined or if input is [], set itest to default
if nargin < 3 
    itest = itestdefault;
elseif isempty(itest)
    itest = itestdefault;
end



% %%%%%%%%%%%%%%%%%%%%%% SET CONTROL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% icut = 5;
% itest = 0; 
% %%%%%%%%%%%%%%%%%%%%%%% SET CONTROL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %use iclean = 0 when no only minimal cleaning of very small clusters is necessary (SEE line 163 below)!!!!<<<<<<<<<<(NORMAL SITUATION)
% %use iclean = 1 when some extra cleaning is necessary 
% %use iclean = 2 when some more extra cleaning is necessary
% %use iclean = -1 when even minimal cleaning is not necessary
% 
% %%% use icut > 0 with Dafi's movie's to remove first few (usually 5) rows (usually contain a lot of noise)
% 
% % use itest = 1 to test the goodness of the cluster detection for the first frame
% % use itest = 0 to make a regular detection run for the clusters in all frames

fprintf(' icut = %d   {icut > 0 for Dafi"s movies to remove first "icut" rows} \n', icut)
fprintf(' iclean = %d  {iclean = 0 is NORMALLY used} \n', iclean)
fprintf(' itest = %d {use itest = 0 for a regular run over a whole movie} \n', itest)

kk=0; %kk is necessary (used as a counter below)
%

% Step 1 : Load all the images 
defaultFileName = []; %HRJ 'tif'

[oriImageName, oriImagePath] = uigetfile('.tif','Please select the first original image'); %HRJ'.tif'
oriImageName = strcat(oriImagePath, oriImageName);
oriImageStackList = getFileStackNames(oriImageName);
total_frame_num = length(oriImageStackList);

%check whether the \maxdata28\ subdirectory exists or not 
%((This is where the data about the clusters is stored))
[success, msg, msgID] = mkdir(oriImagePath,'maxdata283');

if (success ~= 1)
    error(msgID, msg); 
elseif (~isempty(msg))
    fprintf('Subdirectory "maxdata283" already exists.\n');
else
    fprintf('Subdirectory "maxdata283" has been created.\n');
end

%check whether the \images28\ subdirectory exists or not
%((This is where the binary images of the clusters are stored))
[success, msg, msgID] = mkdir(oriImagePath,'images283');

if (success ~= 1)
    error(msgID, msg); 
elseif (~isempty(msg))
    fprintf('Subdirectory "images283" already exists.\n');
else
    fprintf('Subdirectory "images283" has been created.\n');
end


% MovieName = 'Movie283';
% MakeQTMovie ('start', [oriImagePath, filesep, 'maxdata283', filesep, MovieName '.mov']);

ixmax = total_frame_num;
if itest == 1
    ixmax = 1;
end

for ix = 1 : 1 : ixmax
    
    tempname = char(oriImageStackList(ix));
    r = imread(tempname);
    [sy,sx]=size(r);
    %if icut==1
   
    r=r(1+icut:sy,1:sx);
    [sy,sx]=size(r);
    %end
    r=double(r);
    mxo=max(r(:));
    mno=min(r(:));
    r=r*450/mxo; % to renormalize the maximum intensity to a reasonable value.
    %avo=mean(r(:));
    %sdo=std(r(:));
    
    
    k=4;  %We'll calculate details up to order k
    
    rw=wstarck2(r,k,0); % CALCULATES reconstructed image up to order k without kth detail(irecon=0)
    er=r-rw;  % This is the error image
    
    % Now do iterative filtering
    delta=10;
    n=0;
    sig1=std(er(:));
    while delta > 0.002
        
        n=n+1;
        er2=wstarck2(er,k,0);
        rw=rw+er2;
        
        er=r-rw;
        sig2=std(er(:));
        delta=abs(sig2-sig1)/sig2;
        sig1=sig2;
    end
    %     mxerr=max(er(:))
    %     mnerr=min(er(:))
    %     averr=mean(er(:))
    wm=wstarck222(rw,k,1); % calculate the multiscale product of the details
    
    ttmn=min(rw(:));
    rw=rw-ttmn;
    mxw=max(rw(:));
    mnw=min(rw(:));
    
    bs=4; % size of mask for local average= (2*bs+1)by(2*bs+1)
    
    [av,sg]=wlocav(rw,bs); % calculate local average (av) and local standard deviation (sg)
    
    ttav=mean(rw(:));
    ttsg=std(rw(:));
    
    avt=ttav; %+0.5*ttsg;
    
    %create binary image
    
    ia=0; % ia indicates unoccupied (or dark) sites,
    
    ib=1-ia; % ib indicates fluorescent sites
    for i=1:sy
        for j=1:sx
            sd3=av(i,j)+0.5*sg(i,j);
            if ((rw(i,j) >= sd3)&&(rw(i,j)*wm(i,j) >=avt)) % local and global cutoff criteria
                rwb(i,j)=ib;
            else
                rwb(i,j)=ia;
            end
        end
    end
    
    
    rwb=bwmorph(rwb,'clean'); %to get rid of isolated pixels (noise)
    rwb=bwmorph(rwb,'fill'); %to fill up empty internal pixels
    rwb=bwmorph(rwb,'thicken'); %to make larger clusters because of the harsh cutoff criteria above
    if iclean >=0
        rwb=bwmorph(rwb,'spur'); % to remove single pixels 8-attached to clusters
        rwb=bwmorph(rwb,'spur');
        rwb=bwmorph(rwb,'clean');% to remove any remaining isolated pixels after applying spur twice
    end
    if iclean > 0 %Value of iclean is set in line 56 above
        rwb=bwmorph(rwb,'erode');%extra cleaning of small spots
        
        if iclean >=2 
            rwb=bwmorph(rwb,'spur');%extra cleaning of small spots % for extraextra cleaning 
        end
        if iclean >=1 
            
            rwb=bwmorph(rwb,'clean');%extra cleaning of small spots
        end
        rwb=bwmorph(rwb,'thicken');%extra cleaning of small spots
    end
    %     rwb=bwmorph(rwb,'spur');
    %     rwb=bwmorph(rwb,'spur');
    [Lbr,num] = bwlabel(rwb,8); %Label the clusters. Better to use label 8 than label 4
    %fprintf([' num = ',num2str(num),'\n'])
    
    rw=rw*(mxo-mno)/mxw; % normalize maximum intensity to original (maximum intensity-background intensity).
    
    rw=rw.*rwb; %Outside the clusters the intensity is set to zero exactly
    
    ttav=mean(rw(:));
    ttsg=std(rw(:));
    tta1=ttav; %+1.0*ttsg
    ttmx1=max(rw(:));
    ttmn1=min(rw(:));
    
    lm=locmax2d(rw,[9 9]); % find the location of the maxima of the clusters
    [ymx xmx]=find(lm>0); %tta1);
    
    nmax=length(xmx);  % calculate  number of maxima determined by locmax2d
    pnum=0; %initialize number of primary maxima
    snum=0; %initialize number of secondary maxima (if there are two or more maxima in a cluster found by locmax2d)
    lxm=zeros(num,1);  %initialize the number of maxima in each cluster
    intot=zeros(num,1);%initialize the total intensity of each cluster
    csize=zeros(num,1);%initialize the size of each cluster
    xmax=zeros(num,1);% initialize the coordinates of the maxima
    ymax=zeros(num,1);
    yav=zeros(num,1);%initialize the coordinates of the center of intensity
    xav=zeros(num,1);
    labl=zeros(num,1); %NEWHRJ
    for i=1:num
        labl(i)=i;%NEWHRJ
        [yi,xi] = find(Lbr==i);
        kk=kk+1;
        lni=length(yi);
        nn(kk)=lni;
        [ym,xm]= find((Lbr==i)&(lm>0));
        
        lxm(i)=length(xm);
        pnum=pnum+1;
        if lxm(i)==1
            
            ymax(pnum)=ym;
            xmax(pnum)=xm;
            for j=1:lni
                inloc=rw(yi(j),xi(j));
                yav(i)=yav(i)+yi(j)*inloc;
                xav(i)=xav(i)+xi(j)*inloc;
                intot(i)=intot(i)+inloc;
                csize(i)=csize(i)+1;
            end
            xav(i)=xav(i)/intot(i);
            yav(i)=yav(i)/intot(i);
        elseif lxm(i)==0 %if the cluster's maximum has not been located by locmax, search for it locally
            nmax=nmax+1;
            
            tt=0;
            lxm(i)=1;
            for j=1:lni
                inloc=rw(yi(j),xi(j));
                yav(i)=yav(i)+yi(j)*inloc;
                xav(i)=xav(i)+xi(j)*inloc;
                intot(i)=intot(i)+inloc;
                csize(i)=csize(i)+1;
                if inloc>tt
                    tt=inloc;
                    ymx(nmax)=yi(j);
                    xmx(nmax)=xi(j);
                    ymax(pnum)=yi(j);
                    xmax(pnum)=xi(j);
                end
                %fprintf([' i = ',num2str(i),' j=  ',num2str(j),'\n'])
                %ymax(pnum)=ymx(nmax);
                %xmax(pnum)=xmx(nmax);
            end
            xav(i)=xav(i)/intot(i);
            yav(i)=yav(i)/intot(i);
        else  % Take care of the case when there are more than one maximum in each cluster
            
            for j=1:lni
                inloc=rw(yi(j),xi(j));
                yav(i)=yav(i)+yi(j)*inloc;
                xav(i)=xav(i)+xi(j)*inloc;
                intot(i)=intot(i)+inloc;
                csize(i)=csize(i)+1;
            end
            xav(i)=xav(i)/intot(i);
            yav(i)=yav(i)/intot(i);
            tt1=rw(ym(1),xm(1));
            ymax(pnum)=ym(1);
            xmax(pnum)=xm(1);
            for jmx=2:lxm(i)
                tt2=rw(ym(jmx),xm(jmx));
                snum=snum+1;
                if tt2>tt1
                    ymax2(snum)=ymax(pnum); %The weaker maxima  are designated secondary maxima
                    xmax2(snum)=xmax(pnum);
                    ymax(pnum)=ym(jmx); %The maximum with the highest intensity in a cluster is designated the primary max
                    xmax(pnum)=xm(jmx);
                    lxm2(snum)=lxm(i);
                    intot2(snum)=intot(i);
                    csize2(snum)=csize(i);
                    labl2(snum)=i;%NEWHRJ
                    tt1=tt2;
                else
                    ymax2(snum)=ym(jmx);
                    xmax2(snum)=xm(jmx);
                    lxm2(snum)=lxm(i);
                    intot2(snum)=intot(i);
                    csize2(snum)=csize(i);
                    labl2(snum)=i;%NEWHRJ
                end
            end
            
        end
        
    end
    % Add secondary maxima to the end of the primary maxima list
    if snum > 0
        for  imx=1:snum
            pnum=pnum+1;
            ymax(pnum)=ymax2(imx);
            xmax(pnum)=xmax2(imx);
            lxm(pnum)=lxm2(imx);
            intot(pnum)=intot2(imx);
            csize(pnum)=csize2(imx);
            labl(pnum)=labl2(imx);%NEWHRJ
        end
    end
    nmax=length(xmax);
    
    fprintf(['ix =   ',num2str(ix),'   nmax =   ',num2str(nmax),'  num =  ',num2str(num),'\n'])
    
    %rn=ttmx1-rw
    
    inn=zeros(nmax,1);
    
    for ii=1:nmax
        
        inn(ii)=rw(ymax(ii),xmax(ii));
        %fprintf([num2str(ii),' y= ',num2str(ymax(ii)),' x=  ',num2str(xmax(ii)), ' inmax= ',num2str(inn(ii)),' intot=  ',num2str(intot(ii)),' size = ',num2str(csize(ii)),' #max = ',num2str(lxm(ii)),' Label = ',num2str(labl(ii)),'\n'])
        %fprintf([num2str(ii),' yav= ',num2str(yav(ii)),' xav=  ',num2str(xav(ii)), ' #max = ',num2str(lxm(ii)),'\n'])
    end
    if itest==1
        figure(2);
        set(2,'Position',[10 10 1.6*sx 1.6*sy]);
        imshow(r, [], 'notruesize');
        axis on
        hold on;
        %plot(xmax,ymax,'r+');
        plot(xav,yav,'r+');
    end
    %     
    %     MakeQTMovie ('addaxes', gca);
    %     close;
    %     
    name2 = [oriImagePath, filesep, 'maxdata283', filesep,'cdata', num2str(3000+ix),'.mat'];
    
    save (name2,'ymax','xmax','inn','yav','xav','intot','csize','lxm','labl','num','nmax');
    
    % newname = [oriImagePath, filesep, 'mask283', filesep,'cmask', num2str(3000+ix),'.jpg'];
    newname = [oriImagePath, 'images283', filesep,'cmask', num2str(3000+ix),'.jpg'];
    
    % eliminate the frequent error messages by rescaling to uint8
    rw8 = uint8(round(255*(rw/max(rw(:)))));
    
    imwrite(rw8,newname,'jpg','Bitdepth',8,'Mode','lossless','Quality',99);
    %
end

%MakeQTMovie ('finish');

name3 = [oriImagePath, filesep, 'maxdata283', filesep,'nn.mat'];
save(name3,'nn');
figure, hist(nn, 150);
% % 
% % Let the user know we have finished
% name3 = [oriImagePath, filesep, 'maxdata283', filesep,'error.mat'];
% save(name3,'er');
fprintf (1, '\nmain283 finished \n');
