function [idlist,intList]=spotID(slist,opt,dataProperties,projName)
%SPOTID finds corresponding points of subsequent time points in the spotlist (if left away, program asks for it (recommended))
%
% SYNOPSIS  [idlist,intList]=spotID(slist,opt,dataProperties,projName)
%
% INPUT slist : spot list generated with autospfind
%       opt   : optional options structure
%               opt.weight: weight of intensity vs. displacement (0...{1/2}...1)
%               opt.verbose: whether graphs should be displayed ({0}/1/2)
%                   2 only displays intensity graphs
%               opt.save: save or not (0/{1})
%               opt.checkIntensity: (optional) perform intensity check in spotID_findStats (0/{1})
%
% OUTPUT idlist : list of same size as slist and #of elements at time point t same as # of spots in slist(t)
%                 idlist.spot.amp
%                 idlist.spot.coord
%                 idlist.spot.linkup
%                 idlist.spot.linkdown
%                 idlist.spot.color
%                 idlist.spot.background
%                 idlist.linklist: for each time t: 
%                       time/#of spot/color of spot/color of marker/flag/linkup/linkdown/intensity
%
%
% DEPENDS ON: spotID_mappoints, spotID_findStats, bsum2bvec, visualize data with
% labelgui
%   
% c: 18/09/02	Jonas

%----------------initialization


%CONST DEFINITIONS
if nargin<3 | isempty(dataProperties)
    error('no dataProperties found');
%     load_sp_constants;
%     global PIXELSIZE_XY PIXELSIZE_Z  MAXSPOTS;
%     timeLapse=1.001;
else
    
    PIXELSIZE_XY=dataProperties.PIXELSIZE_XY;
    PIXELSIZE_Z=dataProperties.PIXELSIZE_Z;
    MAXSPOTS=dataProperties.MAXSPOTS; 
    timeLapse=dataProperties.timeLapse;
end

if nargin<4|isempty(projName)
    projName='';
end

%get input
if nargin==0|isempty(slist)
    [slistname, pathname] = uigetfile('*.*', 'Select a spotlist (slist)');
    eval(['cd ',pathname,';']); %this command does not work if there is a space in the folder or file name!
    load(slistname);
    %set opt
    verbose=0;
    weight=1/2; 
    saveFile=1;
    findGoodTime=1;
    
elseif nargin==2|nargin==3|nargin==4
    if isfield(opt,'verbose')
        verbose=opt.verbose;
    else
        verbose=0;
    end
    if isfield(opt,'weight')
        weight=opt.weight;
    else
        weight=0.5;
    end
    if isfield(opt,'save')
        saveFile=opt.save;
    else
        saveFile=1;
    end
    if isfield(opt,'checkIntensity')
        findGoodTime=opt.checkIntensity;
    else
        findGoodTime=1;
    end

elseif nargin==1
    verbose=0;
    weight=1/2; %relative weight of intensiti vs. distance for mapping (weight distance/weight intensity; w=[0...1])
    saveFile=1;
    findGoodTime=1;
    
else
        error('more than three or less than zero input arguments')
end

%test if stoopid me has forgotten something
if ischar(weight)
    weight=str2double(weight);
end
%---------------------------------------------------------------------------------------------------


%pixel 2 micron PROGRAM DOES NOT WORK WITH MORE THAN 10 SPOTS!
p2m=ones(10,1)*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];


sigma.xyz=[];
sigma.amp=[];
dist=[];


%read slist into new structure: spots(t).xyz (n*3) coordinates, spots(t).amp (n*1) amplitudes
tmax=size(slist,2);
weight=weight.*ones(tmax,1); %for recalcIdlist
nspots=zeros(tmax,1);
spots(1:tmax)=struct('amp',[],'xyz',[]);
if isfield(slist(1),'sp')
    doFindStats=1;
    for t=1:tmax
        if ~isempty (slist(t).sp)
            spots(t).amp=cat(1,slist(t).sp.amp); %list of amplitudes
            nspots(t)=length(spots(t).amp);  %number of spots in frame t
            if ~isempty(slist(t).statistics)
                spots(t).detectQ=slist(t).statistics.Q;
                spots(t).noise=slist(t).statistics.chi;
            else
                spots(t).detectQ=[];
                spots(t).noise=[];
            end
            spots(t).trackQ=[];
            if nspots(t)>10
                disp(['Warning: more than 10 spots found by detector in frame ',num2str(t),' - can not be handled by spotID']);
                nspots(t)=0;
                spots(t).amp=[];
            else
                spots(t).xyz=cat(1,slist(t).sp.cord).*p2m(1:nspots(t),:); %list of coordinates [um]
            end
        end
    end
else %input came from recalcIdlist
    spots=slist;
    nspots=slist(1).nspots;
    doFindStats=0;
    deltamp=slist(1).deltAmp;
    sigma=slist(1).sigma;
    goodTime=slist(1).goodTime;
    CoMList=slist(1).CoMList;
end

%test if Dom accidentially cut slist so that first frame is a rejected one
tStart=1;
if isempty(spots(1).amp)
    %lookfor first good frame
    done=0;
    while ~done
        tStart=tStart+1;
        if ~isempty(spots(tStart).amp);
            done=1;
        end
    end
end

%-------------------get parameters for analysis--------------------------------
if doFindStats
    [sigma.xyz, sigma.amp, deltamp,goodTime,CoMList,plotData]=...
        spotID_findStats(spots,verbose,findGoodTime,saveFile,timeLapse,projName);
end
%-------------------------------------------------------------------------------

%test if first frame(s) have been discarded by findstats
if goodTime(tStart)==0
    %lookfor first good frame
    done=0;
    while ~done
        tStart=tStart+1;
        if goodTime(tStart)~=0;
            done=1;
        end
    end
end

%ask for weight
if verbose==1
    otherWeight = inputdlg('Do you want to change the weight of intensity to distance (0...1)? Values less than .5 favor intensity Values greater than .5 favor displacement','set weight',1,cellstr(num2str(weight)));
    if ~isempty(otherWeight)
        weight=str2num(otherWeight{1});
        if weight<0 | weight>1
            error('wrong weight (read the instructions)');
        end
    end
end

%----------------calculate which spots link to which

%initialize id.linklist 
%id(t).linklist=[# of frame, # of spot, color(mix), color(pure), flag, linkdown, linkup]
id(1:tmax)=struct('linklist',[]);
id(tStart).linklist=[tStart*ones(nspots(tStart),1),[1:nspots(tStart)]',2.^[0:nspots(tStart)-1]',2.^[0:nspots(tStart)-1]',zeros(nspots(tStart),4)];
%initialize intensitylist: intensities of purecolor-spot [1...n]
intList=zeros(tStart,length(spots(tStart).amp));
intList(tStart,:)=[spots(tStart).amp'];

%launch waitbar
waitbarHandle=mywaitbar(0,[],tmax,'linking tags...');

t1=tStart; %number of first of two frames to compare
nspotsPure(t1)=size(id(t1).linklist,1);
try
    for t2=tStart+1:size(spots,2) %number of second of two frames to compare
        %%%%%%%%%---No mapping if nothing detected at all/if detector did not like the result/if there are more than MAXSPOT-1 spots/if intensity is ok
        %change here to add flags
        if ~isempty(spots(t2).amp)&(spots(t2).amp~=0)&(nspots(t2)<MAXSPOTS)&goodTime(t2)
            
            spots2map=[];
            
            %prepare new step of linklist
            if nspotsPure(t1)<nspots(t2) %new colors are added
                id(t2).linklist=[t2*ones(nspots(t2),1),[1:nspots(t2)]',zeros(nspots(t2),6)];
                newColor=1;
                nspotsPure(t2)=nspots(t2);
            else %there are already enough colors 
                id(t2).linklist=[t2*ones(nspotsPure(t1),1),zeros(nspotsPure(t1),7)];
                newColor=0;
                nspotsPure(t2)=nspotsPure(t1);
            end
            
            %calculate drift correction
            deltaCoM=CoMList(t2,:)-CoMList(t1,:);
            
            %set spots to map to
            for i=1:nspots(t2)
                spots2map(2).amp(i)=spots(t2).amp(i);
                spots2map(2).xyz(i,:)=spots(t2).xyz(i,:)-deltaCoM;
            end
            
            %spots2mapfrom
            for i=1:nspotsPure(t1)
                
                spots2map(1).amp(i)=intList(t1,1+log2(id(t1).linklist(i,4)));
                spots2map(1).xyz(i,:)=spots(t1).xyz(id(t1).linklist(i,2),:);
            end
            
            
            
            nspots2map=[nspotsPure(t1);nspots(t2)];            
            
            %find shortest distance between two points
            [map]=spotID_mappoints(spots2map,nspots2map,mean([sigma.amp(t1),sigma.amp(t2)]),sigma.xyz,sum(deltamp(t1:t2-1)),weight(t2));
            
            dist(t1)=map.dist;
            
            %-------------write results (links&color) into structure
            sameListCell=[]; %init list of spots with same color (for separation)
            switch map.swap
                case 0 %nspots2map(1)>=nspots2map(2)
                    id(t1).linklist(:,7)=map.pL;
                    [c2,up]=sort(map.pL);
                    id(t2).linklist(:,[2,4,6])=[c2',id(t1).linklist(up,4),id(t1).linklist(up,2)];
                    
                    for i=1:nspots(t2) %calculate mix Colors
                        sameIdx=find(id(t2).linklist(:,2)==i);
                        sameCol=sum(id(t2).linklist(sameIdx,4));
                        id(t2).linklist(sameIdx,3)=sameCol*ones(length(sameIdx),1);
                    end
                    
                case 1 %nspots2map(1)<nspots2map(2): a new color has to appear
                    id(t2).linklist(:,[3,4,6])=[id(t1).linklist(map.pL,4),id(t1).linklist(map.pL,4),id(t1).linklist(map.pL,2)];
                    
                    for i=1:nspotsPure(t1)
                        dwn=find(map.pL==i);
                        if length(dwn)==1
                            id(t1).linklist(i,7)=id(t2).linklist(dwn,2);
                        else
                            id(t1).linklist(i,7)=0;
                            sameListCell{end+1}=dwn; %list of the spots in id(t2) with the same color
                        end
                    end
                    
            end
            
            %---check results and improve if needed------------------------------------------------------
            switch sign(nspots(t1)-nspots(t2))
                case 0 %equal number of spots 
                    intList=recalcInt(intList,spots(t2).amp,id(t2).linklist,t1,t2);
                case 1 %more spots in t1: fusion
                    intList=recalcInt(intList,spots(t2).amp,id(t2).linklist,t1,t2);
                    id(t2).linklist(:,5)=1;
                case -1 %less spots in t1: separation 
                    if nspots2map(1)>=nspots(t2) 
                        %no new color
                        intList=recalcInt(intList,spots(t2).amp,id(t2).linklist,t1,t2);
                        %set flag 2
                        id(t2).linklist(:,5)=2;
                        
                    else %add new color
                        %find n spots with same color
                        %loop through sameList: there could be several separations at the same time!
                        for k=1:length(sameListCell)
                            sameList=sameListCell{k};
                            %prepare list of amps and cols for intList&llist-adjustment
                            ampi([1:length(sameList)])=spots2map(2).amp(sameList);
                            coli(1)=id(t2).linklist(sameList(1),4);
                            coli([2:length(sameList)])=max(id(t2).linklist(:,4))*2.^[1:length(sameList)-1]'; 
                            
                            if log2(coli(1))~=round(log2(coli(1))) %test color
                                error('bad mapping (impure spot separation error)');
                            end
                            
                            %calculate intensity ratios
                            ratioi=ampi./sum(ampi);
                            %recalc intensity list with added colors
                            intList(:,log2(coli)+1)=intList(:,log2(coli(1))+1)*ratioi;
                            %recalc linklist
                            for ti=1:t1
                                if ~isempty(id(ti).linklist)
                                    rowIdx=find(id(ti).linklist(:,4)==coli(1)); %find rows with same color
                                    newRow=id(ti).linklist(rowIdx,:);
                                    %adjust mixColor
                                    newRow(3)=newRow(3)+sum(coli(2:end));
                                    %adjust mixColor for all spots
                                    id(ti).linklist(find(id(ti).linklist(:,2)==id(ti).linklist(rowIdx,2)),3)=newRow(3);
                                    %there could be more than one additional colors!
                                    newRow=ones(length(coli)-1,1)*newRow;
                                    for i=1:length(coli)-1
                                        newRow(i,4)=coli(i+1);
                                    end
                                    %write new rows into list
                                    switch (rowIdx==1)+2*(rowIdx==size(id(ti).linklist,1))
                                        case 0
                                            id(ti).linklist=[id(ti).linklist(1:rowIdx,:);newRow;id(ti).linklist(rowIdx+1:end,:)];
                                        case 1
                                            id(ti).linklist=[id(ti).linklist(rowIdx,:);newRow;id(ti).linklist(rowIdx+1:end,:)];
                                        case 2
                                            id(ti).linklist=[id(ti).linklist(1:end,:);newRow];
                                        case 3 %if only one row int ti
                                            id(ti).linklist=[id(ti).linklist(1:end,:);newRow];
                                    end
                                end
                            end
                            
                            %add colors to llist2
                            id(t2).linklist(sameList,4)=coli';
                            id(t2).linklist(sameList,3)=coli';
                            %update intList
                            intList=recalcInt(intList,spots(t2).amp,id(t2).linklist,t1,t2);
                            %set flag 3
                            id(t2).linklist(:,5)=3;
                        end %for-loop
                        
                        %write linkdown for id(t1): if both were sorted according to tag color
                        %sortLL1(:,7)=sortLL2(:,2) -> use sorted indices!
                        [dummy,sortedIdx1]=sortrows(id(t1).linklist,4); %sort LL1
                        [dummy,sortedIdx2]=sortrows(id(t2).linklist,4); %sort LL2
                        id(t1).linklist(sortedIdx1,7)=id(t2).linklist(sortedIdx2,2);
                        
                    end %if-clause
            end
            
            %writing results into structure has to wait till the end - linklist might have to be changed
            t1=t2;
        else
            %timepoint has been rejected. Make sure that this does not result in
            %a plottin error later
            if goodTime(t2)==1
                goodTime(t2) = 0;
            end
        end %if ~isempty(spots(t2).amp)&(spots(t2).amp~=0)&(nspots(t2)<MAXSPOTS)&goodTime(t2)
        mywaitbar(t2/tmax,waitbarHandle,tmax);
    end %for t2=tStart+1:size(spots,2)
catch
    if findstr(lasterr,['Error using ==> get',char(10),'Invalid handle'])
        error('evaluation aborted by user');
    else
        rethrow(lasterror) %or any other action to be executed on an error within the loop
    end
end
close(waitbarHandle);

%---------------write linklist into output structure-----------------------
%prepare for drift-corrected coordinates
coordNoDrift=[];
%loop frames
idlist(1:tmax)=struct('linklist',[],'centroid',[]);
for t=1:size(intList,1)%spots can be longer than intList, if the last frames have not been good
    if ~isempty(id(t).linklist)
        %write linklist
        idlist(t).linklist=id(t).linklist;
        %write intensity into idlist.linklist according to true colors in col4 of linklist
        idlist(t).linklist(:,8)=intList(t,log2(idlist(t).linklist(:,4))+1)';
        %store spot coordinates in linklist (i.e. set tag coordinates)
        idlist(t).linklist(:,9:11)=spots(t).xyz(idlist(t).linklist(:,2),:);
        %store centroid coordinates
        idlist(t).centroid=mean(idlist(t).linklist(:,9:11),1);
        %add chi^2 as 12th col into linklist, but only if not empty (backward-compatibility)
            if ~isempty(spots(t).noise)
                idlist(t).linklist(:,12)=spots(t).noise(idlist(t).linklist(:,2));
            %nse = [nse,spots(t).noise(sp)];
        end
        %sort idlist according to tag color
        idlist(t).linklist = sortrows(idlist(t).linklist,4);
        
        %store Q-matrix from gaussian fit (make sure there is an entry for each tag!)
        detQ = [];
        traQ = [];
        nse  = [];
        %QMatrix will from now on be sorted according to tag color, no
        %longer according to spot number!
        for i=1:size(idlist(t).linklist,1)
            sp = idlist(t).linklist(i,2);
            detQ = blkdiag(detQ,spots(t).detectQ( (sp-1)*3+1:sp*3,(sp-1)*3+1:sp*3 ) );
            
        end
        if ~isempty(spots(t).trackQ) & (size(spots(t).trackQ) ~= 3*size(idlist(t).linklist,1))
            for i=1:size(idlist(t).linklist,1)
                sp = idlist(t).linklist(i,2);
                traQ = blkdiag(traQ,spots(t).trackQ( (sp-1)*3+1:sp*3,(sp-1)*3+1:sp*3 ) );
            end
        else
            traQ = spots(t).trackQ;
        end
        idlist(t).info.detectQ_Pix=detQ;
        %idlist(t).info.noise=nse;
        idlist(t).info.trackQ_Pix=traQ;
        
        %do not forget the trackerMessage!
        if isfield(slist,'trackerMessage')
            idlist(t).info.trackerMessage = slist(t).trackerMessage;
        end
    end
end




%----------------------write data statistics

%write explanation for linklist
idlist(1).stats.help{1,1}='columns of linklist';
idlist(1).stats.help{2,1}='1: timespot #';
idlist(1).stats.help{3,1}='2: spot #';
idlist(1).stats.help{4,1}='3: spot color';
idlist(1).stats.help{5,1}='4: tag color';
idlist(1).stats.help{6,1}='5: flag';
idlist(1).stats.help{7,1}='6: linkup to spot #';
idlist(1).stats.help{8,1}='7: linkdown to spot #';
idlist(1).stats.help{9,1}='8: intensity';
idlist(1).stats.help{10,1}='9-11: x/y/z-coordinates in um (Image Coordinates!)';
idlist(1).stats.help{11,1}='12: chi^2 of the spot on which the tag is located';

%write labellist
idlist(1).stats.labellist{1,1}='spb1';
idlist(1).stats.labellist{2,1}='cen1';
idlist(1).stats.labellist{3,1}='spb2';
idlist(1).stats.labellist{4,1}='cen2';

%write color<->label
idlist(1).stats.labelcolor(1:size(idlist(tStart).linklist,1),1)=cellstr('?');

%deltaAmplitude from findStats
idlist(1).stats.deltAmp=deltamp;

%sigmas for displacement and amplitudes
idlist(1).stats.sigma=sigma;

%maxColor: max value of color
idlist(1).stats.maxColor=sum(idlist(tStart).linklist(:,4))+1;

%weight for spotID_mappoints
idlist(1).stats.weight=weight;

%date of creation for version check
idlist(1).stats.created=date;

%write idlist status
idlist(1).stats.status{1}=[date, ': pure&original idlist'];


%------------------rerun spotID if not recalc
% if doFindStats==1
%     idlist=recalcIdlist(idlist,tStart,[],dataProperties);
%     %rebuild intList
%     for t=1:tmax
%         if ~isempty(idlist(t).linklist)
%             intList(t,log2(idlist(t).linklist(:,4))+1)=idlist(t).linklist(:,8)';
%         end
%     end
% end

%----------------------if verbose: plot intList
if (verbose~=0|saveFile==1) & (findGoodTime == 1)
    gT=find(goodTime);
    if timeLapse==1.001 %no data properties, time not known
        xData=gT;
    else %time is known - still show only timepoints and not absolute time: better for checking.
        %xData=(gT-1)*timeLapse;
        xData=gT;
    end
    figure(plotData.figH);
    hold on;
    
    %draw in "true" colors
    cMap=hsv(64); %init colorMap
    cMapFact=size(cMap,1)/idlist(1).stats.maxColor;
    labelColor=cMap(bsum2bvec(idlist(1).stats.maxColor-1)*cMapFact,:);
    set(gca,'ColorOrder',labelColor);
    
    %draw the line
    spotIntLineH = plot(xData,intList(gT,:));
    
    %get figure handles other than the standard figure handles
    figHandles = guidata(plotData.figH);
    figHandles.spotIntLineH = spotIntLineH;
    %store handles of line
    guidata(plotData.figH, figHandles);
    
    if saveFile==1
        saveas(plotData.figH,['intFigure',nowString,'.fig']);
    end
    if verbose==0
        close(plotData.figH);
        %close(plotData.sigH);
    end
end

%save lists if option is set
idlistFileName=['idlist-',nowString];
intListFileName=['intList-',nowString];


if saveFile
    %save file
    save(idlistFileName,'idlist');
    save(intListFileName,'intList');
end

idlist(1).stats.name = dataProperties.name;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function intList=recalcInt(intListSub,amp,linklist,t1,t2)
%add intensities(t2) to intList

for i=1:max(linklist(:,2))
    %read pure color indices
    iLidx=log2(linklist(find(linklist(:,2)==i),4))'+1;
    %get previous intensities
    iL=intListSub(t1,iLidx);
    %calc new intensities
    intListSub(t2,iLidx)=iL*amp(i)/sum(iL);
end

intList=intListSub;

