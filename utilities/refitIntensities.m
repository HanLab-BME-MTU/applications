function [idlist,dataProperties] = refitIntensities(idlist,dataProperties,movieName,fitPsf,fitPos)
%REFITINTENSITIES fits spot intensities in chromdyn movies
%
% SYNOPSIS: [idlist,dataProperties] = refitIntensities(idlist,dataProperties,movieName,fitPsf,fitPos)
%
% INPUT all input is optional. If the first three arguments are empty, code
%       will ask for directory containing idlists to fit (and save them)
%		idlist: idlist (normally after running tracker)
%		dataProperties: see defaultDataProperties
%		movieName: filename with path of movie
%		fitPsf: if 1, code will attempt to recalculate psf-width (default: 0)
%		fitPos: if 1, code will attempt to re-fit positions, too. (default: 0)
%
% OUTPUT idlist: updated idlist if idlist was supplied. Otherwise, list of
%                updated idlists
%		 dataProperties: updated dataProperties (if fitPsf)
%
% REMARKS
%
% created with MATLAB ver.: 7.5.0.342 (R2007b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 08-Dec-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===========================
%% CHECK INPUT
%===========================

% defaults & pre-assignments
loadData = false;
def_fitPsf = false;
def_fitPos = false;

if nargin < 3 || isempty(idlist) || isempty(dataProperties) || isempty(movieName)
    % load data
    loadData = true;
end

% assume idlist etc. are ok.

% check for fit-vars
if nargin < 4 || isempty(fitPsf)
    fitPsf = def_fitPsf;
end
if nargin < 5 || isempty(fitPos)
    fitPos = def_fitPos;
end

%============================


%============================
%% LOAD DATA
%============================

if loadData
    % read idlistList
    idlistList = loadIdlistList(cdBiodata(4),[],[],1);

else
    % put data into idlistList. Fields needed are:
    % - idlist
    % - loadStruct
    % - dataProperties
    [dummy,dummy,loadStruct] = cdLoadMovie({movieName,'corrected'},'',-1);
    idlistList = struct('idlist',idlist,'dataProperties',dataProperties,...
        'loadStruct',loadStruct);
end

% count
nIdlists = length(idlistList);

%============================


%============================
%% FIT IDLIST
%============================

if loadData
    progressText(0,'Fitting intensities')
end

for iIdlist = 1:nIdlists

    out = idlist2slist(idlistList(iIdlist).idlist,idlistList(iIdlist).dataProperties);

    fitStruct = struct('slist',out,...
        'dataProperties',idlistList(iIdlist).dataProperties,...
        'movieDir',idlistList(iIdlist).loadStruct.moviePath,...
        'rawMovieName',idlistList(iIdlist).loadStruct.movieName);

    %spotIdx = catStruct(1,'out.sp.idxL');

    % for isolated spots, position should always be ok (and, in
    % principle, amplitude, though it will depend on sigma). To get
    % correct sigma, fit width of isolated spots first
    if fitPsf
        [dummy,spbIdx] = ...
            ismember({'spb1';'spb2'},idlistList(iIdlist).idlist(1).stats.labelcolor);
        fitStruct.spbIdx = spbIdx;
        fitStruct = cdFitPsf(fitStruct);
    end

    % check type. If we come from idlist_L, we want to adjust
    % positions, too
    if ~fitPos || findstr(idlistList(iIdlist).idlist(1).stats.idname,'track_L')
        % only fit amplitude
        gaussFit = cdFitGauss(fitStruct,{'a','b'});
    else
        % if it was idlisttrack: find idlist_L
        if findstr(idlistList(iIdlist).idlist(1).stats.idname,'track')
            if strcmp(idlistList(iIdlist).type(end),'2')
                load(idlistList(iIdlist).dataFileName,'idlist_L2')
                idlistList(iIdlist).idlist = idlist_L2;
                fitStruct.idlist = idlist_L2;
            else
                load(idlistList(iIdlist).dataFileName,'idlist_L')
                idlistList(iIdlist).idlist = idlist_L;
                fitStruct.idlist = idlist_L;
            end
        end
        % fit amplitude and position
        gaussFit = cdFitGauss(fitStruct,{'x1','x2','x3','a','b'});
    end

    pix2mu = [fitStruct.dataProperties.PIXELSIZE_XY,...
        fitStruct.dataProperties.PIXELSIZE_XY,...
        fitStruct.dataProperties.PIXELSIZE_Z];

    % write gaussFit back into idlist - transform!
    for t=1:length(idlistList(iIdlist).idlist)
        % positions - if not fitted, they should have stayed the same
        if ~isempty(idlistList(iIdlist).idlist(t).linklist)
            spotIdx = idlistList(iIdlist).idlist(t).linklist(:,2);
            idxL = spotIdx>0; % good spotIdx
            n = sum(idxL);
            idlistList(iIdlist).idlist(t).linklist(idxL,9:11) = ...
                gaussFit(t).coords(:,1:3).*repmat(pix2mu,n,1);
            % amplitudes
            idlistList(iIdlist).idlist(t).linklist(idxL,8) = gaussFit(t).coords(:,4);
            % chi2
            idlistList(iIdlist).idlist(t).linklist(idxL,12) = gaussFit(t).chi2(:);
            % add background
            idlistList(iIdlist).idlist(t).linklist(idxL,16) = gaussFit(t).coords(:,end);
            % Q
            if fitPos
                % update Q

            end
        end
    end
    
    if loadData
        progressText(iIdlist/nIdlists)
    end

end % loop idlists

% check whether we need to save
if loadData

    for iIdlist = 1:nIdlists
        % save into data file
        save(idlistList(iIdlist).dataFileName,idlistList(iIdlist).idlist(1).stats.idname,'-append');
        % save outside data file - later
    end

else
    % return outputs
    idlist = idlistList.idlist;
    dataProperties = idlistList.dataProperties;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = idlist2slist(idlist,dataProperties)

% make a new slist from th eidlist
pix2mu = [dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];

out(1:length(idlist)) = struct('sp',[]);

for t=1:length(idlist)
    if ~isempty(idlist(t).linklist)
        % keep the order of spots of the idlist; don't worry about 'true'
        % spot numbers
        spotIdx = find(idlist(t).linklist(:,2)>0);
        % tracker may assign new spot numbers. Therefore initialize with a
        % common bg
        background = mean(idlist(t).linklist(idlist(t).linklist(:,3)==0,12));
        for i=1:length(spotIdx)
            out(t).sp(i).cord = idlist(t).linklist(spotIdx(i),9:11)./pix2mu;
            out(t).sp(i).amp = idlist(t).linklist(spotIdx(i),8);
            out(t).sp(i).bg = background;
            out(t).sp(i).idxL = spotIdx(i); % row of idlist corresponding to the spot
        end
    end
end