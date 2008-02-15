function [goodIdlist,errorMessage,goodTimes] = checkIdlist(idlist,check,askOptions)
%CHECKIDLIST makes sure that the latest idlist is being used
%
% SYNOPSIS: goodIdlist = checkIdlist(idlist,check)
%
% INPUT idlist: any idlist
%		check : (opt) check to perform. Input can be an integer, a string
%                or a cell array. Depending on the kind of check, it can be
%                performed on the overall idlist, or on individual frames.
%                Some tests accept additional parameters. If you do not
%                want to use the defaults, specify them with help of a cell
%                array. It is also possible to choose the checks
%                interactively, by setting check to 'ask'.
%
%                STRING Strings will be evaluated as is with the eval
%                  command. An example of a check-string is
%                  '~any(strcmp(idlist(1).stats.labelcolor,''cen1*''))'
%                INTEGER Integers point to a specific test:
%                    0: Always true
%                    1: (default) check whether any kind of idlist2
%                    2: check whether there is an inconsistency in the
%                       idlist
%                    3: No cen1* in idlist (same as condition above)
%                    4: Only choose idlists with N1:N2 tags found. Default:
%                       [3:4]
%                    5: Only choose frames with N1:N2 tags found. Default:
%                       [4 4]
%                    6: Only choose frames where avg(cenAmp)/avg(spbAmp) is
%                       between R1 and R2. [0 0.75]
%                    7: Only keep frames where TAG1 and TAG2
%                       - have been found
%                       - arent't fused to each other
%                       - aren't fusing to or separating from another tag.
%                         For example, if tag1 fuses to tag3 in frames 2 to
%                         7, frames 2 and 7 will be removed.
%                       Default: {'spb1', 'cen1'} (needs 2 cols in cell
%                                array)
%                    8: Only keep frames where the spindle pole body
%                       separation is between R1 and R2 um. Def.: [1.4 1.7]
%                       (needs 2 cols in cell array)
%                       If there are only two spots in the movie,  the
%                       distance between the two will be considered as
%                       inter-spb distance, irrespective of labels.
%                    9: No '?' in labels.
%                   10: Only keeps frames wehre the ratio of the projection
%                       of spb1-cen1 onto spb1-spb2 to spb1-spb2 is betwen
%                       R1 and R2. Note: This automatically discards
%                       idlists with less than 3 tags found. Def.: [0 0.25]
%                CELL A cell array is used to specify multiple checks
%                  and/or use values other than the defaults. The
%                  dimensions of the cell array are nChecks-by-m, where m
%                  can be 1 to 3. The first column of the cell array lists
%                  the test, either as string or as integer (the column,
%                  but not a cell, may contain both). For integers, the
%                  second column specifies the optional additional
%                  parameters. An example of a cell is:
%                  {'~any(strcmp(idlist(1).stats.labelcolor,''cen1*''))','';...
%                   5, 3:4;...
%                   7, 'spb2', 'cen2'}
%                   An empty cell array as input will always evaluate to
%                   true.
%                'ask' lets the user choose via a GUI which checks to use.
%                  The output of the choice can be made to apply for all
%                  subsequent idlists that are loaded using persistent
%                  variables.
%                  WARNING: if checkIdlist is called with 'ask', the
%                  calling function has to end with
%                    clear checkIdlist
%                  Otherwise, the choice of checks will apply the next time
%                  checkIdlist is called with 'ask' during the same Matlab
%                  session.
%       askOptions : optional structure with options for the dialog window.
%                  .choiceList : list of checks to presented as a choice
%                       (see above for numbering). Default: 3:8
%                  .choiceDefaults : nChoices-by-2 cell array with defaults
%                  .checkCell : Additional checks to be run. Needs to be in
%                       cell-array format. These checks will not be
%                       presented to the user.
%
% OUTPUT goodIdlist: logical 1 or 0, depending on whether the idlist passed
%                       the required check
%        errorMessage : message in case of ~goodIdlist
%        goodTimes : logical vector of the same length as idlist containing
%                    1 for timepoints that pass the check. If all(goodTimes
%                    == 0), goodIdlist will be set to 0. To remove the bad
%                    entries in idlist, use
%                      idlist(~goodTimes).linklist = [];
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 27-Apr-2006
% greatly expanded in January 2008 by Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==================
%% SET DEFAULTS
%==================

% if you add checks, don't forget to update the subfunction

def_nTagsInIdlist = 3:4;
def_nTagsInFrame  = 4;
def_ampThreshold = [0 0.8];
def_goodTags = {'spb1','cen1'};
def_spbSeparation = [1.4 1.7];
def_choiceList = [3:6,8:10];
def_distanceRatio = [0,0.25];

choiceDefaults = cell(10,2);
choiceDefaults{4,1} = def_nTagsInIdlist;
choiceDefaults{5,1} = def_nTagsInFrame;
choiceDefaults{6,1} = def_ampThreshold(1);
choiceDefaults{6,2} = def_ampThreshold(2);
choiceDefaults(7,1:2) = def_goodTags;
choiceDefaults{8,1} = def_spbSeparation(1);
choiceDefaults{8,2} = def_spbSeparation(2);
choiceDefaults{10,1} = def_distanceRatio(1);
choiceDefaults{10,2} = def_distanceRatio(2);

% define persistent_check here, because it cannot be inside an if-statement
 persistent persistent_check

%==================
%% TEST INPUT
%==================

if nargin < 1 || isempty(idlist)
    error('checkIdlist needs at least an idlist as input!')
end

% preassign output. Initialize goodTimes to all true, because the overall
% tests don't care about frames, and have the frame-by-frame tests get a
% tmpGoodTimes that is then combined with the existing one via AND.
% Initialize goodIdlist as true, too, to facilitate combining tests.
% preassign here, so that we can quit if no test requested.
goodTimes = true(length(idlist),1);
goodIdlist = true;
errorMessage = '';



% nargin can't be smaller than 2 with persistent check, because it has to
% be 'ask' initially to allow creating a non-persistent variable
if nargin < 2
    % default check is 1
    check = 1;
end

if isempty(check)
    if iscell(check)
        % empty check. Always true.
        return
    else
        % default check is one
        check = 1;
    end
end

%% ask for input
if ischar(check) && strcmp(check,'ask')

    % check for persistent variables
   
    if ~isempty(persistent_check)
        check = persistent_check;
    else

        % check for askOptions
        if nargin < 3
            askOptions = [];
        end
        if isfield(askOptions,'checkCell')
            checkCell = askOptions.checkCell;
            scc = size(checkCell);
            if scc(2) ~= 3
                checkCellTmp = checkCell;
                checkCell = cell(scc(1),3);
                checkCell(:,1:min(scc(2),3)) = checkCellTmp(:,1:min(scc(2),3));
            end
        else
            % have at least check 0, so that an acceptAll will make the ask
            % window disappear
            checkCell = {0,'',''};
        end
        if isfield(askOptions,'choiceList')
            choiceList = askOptions.choiceList;
        else
            choiceList = def_choiceList;
        end
        % make sure checks in checkCell don't show up in choiceList
        if ~isempty(checkCell)
            c = [checkCell{:,1}];
            badChoice = ismember(choiceList,c);
            if isfield(askOptions,'choiceDefaults')
                askOptions.choiceDefaults(badChoice,:) = [];
            end
            choiceList(badChoice) = [];
        end
        if isfield(askOptions,'choiceDefaults')
            choiceDefaults = askOptions.choiceDefaults;
        else


            choiceDefaults = choiceDefaults(choiceList,:);
        end



        % prepare input for checkIdlistGui
        choiceCell = getChoiceCell;

        % add defaults
        choiceCell(choiceList(:)+1,[4 6]) = choiceDefaults;
        % retain good choices. Don't forget to keep empty choice
        choiceCell = choiceCell([1;choiceList(:)+1],:);

        % submit
        [check, applyToAll] = checkIdlistGui(choiceCell, idlist(1).stats.name);

        if ~iscell(check) && check == -1
            error('checkIdlist aborted')
        end

        % add checkCell
        check = [check;checkCell];

        % make persistent if requested
        if applyToAll
            persistent_check = check;
        end
    end

end

%% continue testing input

% check for cell input
if iscell(check)
    % check for empty cell
    if isempty(check)
        return
    end
    % get a vector of checks. User has to supply correct r/c, because I
    % can't test for it without huge effort.
    checkCell = check;
    check = check(:,1);
    nChecks = length(check);
else
    % turn check into a cell input
    nChecks = 1;
    checkCell = [];
    check = {check};
end



%===================
%% CHECK IDLIST
%===================

% loop checks
iCheck = 0;
while goodIdlist && iCheck < nChecks
    iCheck = iCheck +1;
    switch check{iCheck}

        case 0
            % true
            goodIdlist = goodIdlist && true;

        case 1 % check for idlist2

            % old idlists didn't have a field 'intFit'
            goodIdlist = goodIdlist && ...
                isfield(idlist,'stats') && isfield(idlist(1).stats,'intFit');
            if ~goodIdlist
                errorMessage = 'Old idlist!';
            end


        case 2 % check for inconsistencies

            % Bug #164 - not the same number of tags as in labelcolor

            % find goodTags from linklist
            linklists = cat(3,idlist.linklist);
            % badTagIdx = find(ismember(linklists(:,5,1),[2,3]));
            % goodTagIdx = missingIndices(badTagIdx,size(linklists,1));
            % compare to labelcolor
            if size(linklists,1) ~= length(idlist(1).stats.labelcolor);%length(goodTagIdx) ~= length(idlist(1).stats.labelcolor)
                errorMessage = ...
                    'Inconsistency in idlist detected: wrong number of tag labels!';
                goodIdlist = false;
            end

            for t = squeeze(linklists(1,1,:))'
                if isfield(idlist(t),'trackInit') && ...
                        ~isempty(idlist(t).trackInit) && ...
                        ~all(ismember(idlist(t).trackInit(:,1),idlist(t).linklist(:,4)))
                    errorMessage = 'TrackInit points to non-existent tag!';
                    goodIdlist = false;
                    break
                end
            end


        case 3 % avoid cen1*

            goodIdlist = goodIdlist && ...
                ~any(strcmp(idlist(1).stats.labelcolor,'cen1*'));

            if ~goodIdlist
                errorMessage = 'Idlist contains cen1*';
            end


        case 4 % check for the right number of good tags

            % remove first all single-occurence tags, so that we're only
            % left with good ones. LG_deleteTag needs goodTimes. Don't call
            % it that way, though.
            if checkIdlist(idlist,1)
                gt = catStruct(1,'idlist.linklist(1)');
                idlistTmp = LG_deleteSingleOccurences(idlist,gt);
            else
                idlistTmp = idlist;
            end

            % now read out number of tags
            nTags = length(idlistTmp(1).stats.labelcolor);

            % check whether there is something else than the default
            if isempty(checkCell) || size(checkCell,2) < 2 || isempty(checkCell{iCheck,2})
                nTagsInIdlist = def_nTagsInIdlist;
            else
                nTagsInIdlist = checkCell{iCheck,2};
            end

            % evaluate condition
            goodIdlist = goodIdlist && any(nTagsInIdlist == nTags);

            % write error message
            if ~goodIdlist
                errorMessage = sprintf('Idlist contains %i tags instead of (%s\b)',nTags,sprintf('%i/',nTagsInIdlist));
            end

        case 5 % check for the right number of good tags in each frame

            % remove first all single-occurence tags, so that we're only
            % left with good ones. LG_deleteTag needs goodTimes. Don't call
            % it that way, though.
            if checkIdlist(idlist,1)
                gt = catStruct(1,'idlist.linklist(1)');
                idlistTmp = LG_deleteSingleOccurences(idlist,gt);
            else
                idlistTmp = idlist;
            end


            % find target number of tags
            if isempty(checkCell) || size(checkCell,2) < 2 || isempty(checkCell{iCheck,2})
                nTagsInFrame = def_nTagsInFrame;
            else
                nTagsInFrame = checkCell{iCheck,2};
            end

            % loop through idlist, fill goodTimes with zeros. Empty idlist
            % counts as bad. For safety reasons, make a temporary
            % goodTimes-vector before folding it with the output
            gtTmp = false(size(goodTimes));
            for t = 1:length(idlistTmp)
                % check for empty idlist, check for the proper number of
                % detected/tracked tags
                if ~isempty(idlistTmp(t).linklist) && any(sum(idlistTmp(t).linklist(:,3)~=1)==nTagsInFrame)
                    gtTmp(t) = true;
                end
            end

            % combine gtTmp with goodTimes
            goodTimes = goodTimes & gtTmp;


        case 6 % check for the right amplitudes

            % find the tagIndices
            spbIdx = strmatch('spb',idlist(1).stats.labelcolor);
            cenIdx = strmatch('cen',idlist(1).stats.labelcolor);

            % find the target value
            if isempty(checkCell) || size(checkCell,2) < 2 || isempty(checkCell{iCheck,2})
                ampThreshold = def_ampThreshold;
            else
                ampThreshold = checkCell{iCheck,2};
            end

            % loop through idlist and check
            gtTmp = false(size(goodTimes));
            for t = 1:length(idlist)
                if ~isempty(idlist(t).linklist)
                    if goodTimes(t)
                        % read amps, discard estimated tags, secondary fusions
                        spbAmp = idlist(t).linklist(spbIdx,8);
                        spbAmp(ismember(idlist(t).linklist(spbIdx,3),[1 4])) = [];
                        cenAmp = idlist(t).linklist(cenIdx,8);
                        cenAmp(ismember(idlist(t).linklist(cenIdx,3),[1 4])) = [];

                        % check target
                        if mean(cenAmp)/mean(spbAmp) > ampThreshold(1) && mean(cenAmp)/mean(spbAmp) < ampThreshold(1)
                            gtTmp(t) = true;
                        end
                    end
                end

            end % loop idlist

            % combine gtTmp with goodTimes
            goodTimes = goodTimes & gtTmp;


        case 7 % make sure tags are there and there are no fusion-related jumps

            %             7: Only keep frames where TAG1 and TAG2
            %             - have been found
            %             - arent't fused to each other
            %             - aren't fusing to or separating from another tag.
            %             For example, if tag1 fuses to tag3 in frames 2 to
            %                 7, frames 2 and 7 will be removed.
            %                 Default: {'spb1', 'cen1'} (needs 2 cols in cell
            %                 array)



            % check for tag names
            if isempty(checkCell) || size(checkCell,2) < 3 || isempty(checkCell{iCheck,2})
                goodTags = def_goodTags;
            else
                goodTags = checkCell(iCheck,2:3);
            end

            % find tags
            tag1Idx = strmatch(goodTags{1},idlist(1).stats.labelcolor);
            tag2Idx = strmatch(goodTags{2},idlist(1).stats.labelcolor);
            tagIdx = [tag1Idx,tag2Idx];

            nTags = length(idlist(1).stats.labelcolor);

            % set all goodTimes to false
            gtTmp = false(size(goodTimes));

            % check for new idlist
            newId = checkIdlist(idlist);

            % loop through frames, check
            fused = [0,0];
            for t = 1:length(idlist)

                if goodTimes(t)
                    if ~isempty(idlist(t).linklist)


                        % check for both tags existing
                        if newId
                            if ~any(ismember(idlist(t).linklist(tagIdx,3),[1,4]))
                                gtTmp(t) = true;
                            end
                        else
                            if ~any(idlist(t).linklist(tagIdx,2)==0)
                                gtTmp(t) = true;
                            end
                        end

                        % check for fusion
                        tagCoords = idlist(t).linklist(tagIdx,9:11);

                        if all(tagCoords(1,:) == tagCoords(2,:))
                            % both fused. Discard
                            gtTmp(t) = false;
                        else
                            % check for fusion to other tags

                            % tag1
                            not1 = true(nTags,1);
                            not1(tag1Idx) = false;
                            if any(ismember(tagCoords(1,:),idlist(t).linklist(not1,9:11),'rows'))
                                if fused(1) == 0
                                    % not fused before. Discard
                                    gtTmp(t) = false;
                                    % remember time
                                    fused(1) = t;
                                else
                                    % remember time
                                    fused(1) = t;
                                end
                            else
                                % check for prior fusion
                                if fused(1) == 0
                                    % all is well
                                else
                                    % remove last fusion
                                    gtTmp(fused(1)) = false;
                                    % remember not fused
                                    fused(1) = 0;
                                end
                            end

                            % tag2
                            not2 = true(nTags,1);
                            not2(tag2Idx) = false;
                            if any(ismember(tagCoords(2,:),idlist(t).linklist(not2,9:11),'rows'))
                                if fused(2) == 0
                                    % not fused before. Discard
                                    gtTmp(t) = false;
                                    % remember time
                                    fused(2) = t;
                                else
                                    % remember time
                                    fused(2) = t;
                                end
                            else
                                % check for prior fusion
                                if fused(2) == 0
                                    % all is well
                                else
                                    % remove last fusion
                                    gtTmp(fused(2)) = false;
                                    % remember not fused
                                    fused(2) = 0;
                                end
                            end
                        end % if fused tags

                    end
                end % if goodTimes
            end % loop time

            % combine gtTmp with goodTimes
            goodTimes = goodTimes & gtTmp;


        case 8 % check spb separation

            % find the tagIndices
            spbIdx = strmatch('spb',idlist(1).stats.labelcolor);

            % if there are only two spots, spbIdx is all of them
            if length(idlist(1).stats.labelcolor) == 2
                spbIdx = [1;2];
            end

            % find the target value
            if isempty(checkCell) || size(checkCell,2) < 3 || isempty(checkCell{iCheck,2})
                spbSeparation = def_spbSeparation;
            else
                spbSeparation = [checkCell{iCheck,2:3}];
            end

            % loop through idlist and check
            gtTmp = false(size(goodTimes));
            for t = 1:length(idlist)
                if goodTimes(t)
                    if ~isempty(idlist(t).linklist)
                        % read spb-coord
                        spbCoord = idlist(t).linklist(spbIdx,9:11);
                        % calc distance
                        spbDeltaCoord = norm(diff(spbCoord));

                        % check target
                        if spbDeltaCoord > spbSeparation(1) && spbDeltaCoord < spbSeparation(2)
                            gtTmp(t) = true;
                        end
                    end
                end % if goodTimes

            end % loop idlist

            % combine gtTmp with goodTimes
            goodTimes = goodTimes & gtTmp;

        case 9 % no unlabeled good tags
            % remove first all single-occurence tags, so that we're only
            % left with good ones. LG_deleteTag needs goodTimes. Don't call
            % it that way, though.
            gt = catStruct(1,'idlist.linklist(1)');
            if checkIdlist
                idlistTmp = LG_deleteSingleOccurences(idlist,gt);
            else
                idlistTmp = idlist;
            end

            goodIdlist = goodIdlist && ...
                ~any(strcmp(idlistTmp(1).stats.labelcolor,'?'));

            if ~goodIdlist
                errorMessage = 'Idlist contains ''?''';
            end

        case 10 % ratio of projection of s1c1 onto s1s2

            % check for correct number of spots
            goodIdlist = checkIdlist(idlist,4);

            if goodIdlist

                % check for target value
                if isempty(checkCell) || size(checkCell,2) < 3 || isempty(checkCell{iCheck,2})
                    distanceRatio = def_distanceRatio;
                else
                    distanceRatio = [checkCell{iCheck,2:3}];
                end

                % find tags
                [tagExists,tagOrder] = ismember({'spb1','cen1','cen2','spb2'}, idlist(1).stats.labelcolor);


                % calculate tagOrder for distance calculation
                idTagOrder = tagOrder(tagOrder>0);
                idTagOrder = idTagOrder - min(idTagOrder) + 1;

                if checkIdlist( idlist,1)
                    % new idlist - dataProperties are needed for sigmas,
                    % which is not relevant for now
                    [ distance,  distanceUnitVectors] = ...
                        idlist2distMat( idlist,  defaultDataProperties,[],[],idTagOrder);
                else

                    % calculate distances, distanceVectors without idlist2distMat - we need to
                    % be able to use old idlists. Since we don't worry about uncertainties,
                    % it's still fairly straightforward
                    linklists = cat(3, idlist.linklist);
                    linklists = linklists(idTagOrder,:,:);
                    goodTimes = squeeze(linklists(1,1,:));
                    nTimepoints = length( idlist);
                    nTags = length(idTagOrder);
                    distance = repmat(NaN,[nTags,nTags,nTimepoints]);
                    distanceUnitVectors = repmat(NaN,[nTags,nTags,nTimepoints,3]);
                    for t = goodTimes'
                        % get distance matrix, distanceVectorMatrix. Divide
                        % distanceVectorMatrix by distance to get normed vectors
                        % use the same ordering as idlist2distMat
                        [ distance(:,:,t),distanceVectorMatrix] =...
                            distMat(linklists(:,9:11,linklists(1,1,:)==t));
                        distanceUnitVectors(:,:,t,:) = ...
                            permute(...
                            distanceVectorMatrix./repmat( distance(:,:,t),[1,1,3]),[1,2,4,3]);
                    end
                end

                % store index, time
                distList = [squeeze( distance(2,1,:)),squeeze( distance(3,1,:))];


                % get scalar product of vectors for projection of s1c1 onto
                % s1s2
                s1c1vec = squeeze( distanceUnitVectors(2,1,:,:));
                s1s2vec = squeeze( distanceUnitVectors(3,1,:,:));
                proj = sum(s1c1vec.*s1s2vec,2);


                distList(:,3) = ( distList(:,1).*proj)./ distList(:,2);
                if sum( distList(:,3)>0.5) > 0.5 * sum(~isnan(proj))
                    distList(:,3) = 1- distList(:,3);
                end

                % update goodTimes
                goodTimes = goodTimes & (distList(:,3)>distanceRatio(1) & distList(:,3)<distanceRatio(2));

            end

        otherwise
            % check for string
            if ischar(check{iCheck})
                goodIdlist = goodIdlist && ...
                    eval(check{iCheck});
            else
                error('check %i has not been defined yet',check)
            end
    end % switch check{iCheck}

    % check whether there is any goodTime
    if ~any(goodTimes)
        goodIdlist = false;
        errorMessage = 'No good frame detected';
    end

end % while goodIdlist && iCheck < nChecks




function choiceCell = getChoiceCell
% subfunction listing the possible choices

% name (in pd), choice number, inputType, default, title, default, title
% inputType: 0 - none
%            1 - numeric array
%            2 - string
%           11 - 2 numeric arrays
%           12 - numeric/string
%           21 - string/numeric
%           22 - 2 strings
choiceCell = {'Please choose',0, 0, '', '', '', '';...
    'is idlist2?' ,1, 0, '', '', '', '';...
    'is inconsistent?',2, 0, '', '', '', '';...
    'no cen1*',3, 0, '', '', '', '';...
    'idlist with N1:N2 tags',4, 1, '', 'N1:N2', '', '';...
    'frames with N1:N2 spots',5, 1, '', 'N1:N2', '', '';...
    'frames w/ avg(cenAmp)/avg(spbAmp)<R',6, 11, '', 'maxRatio', '', '';...
    'frames with tag1,tag2',7, 22, '', 'tag1', '', 'tag2';...
    'frames with minD<d(spb)<maxD (um)',8, 11, '', 'minD (um)', '', 'maxD (um)';...
    'no unlabeled good tags (''?'')',9,0,'','','','';...
    'frames w/ minR<s1c1/s1s2<maxR',10,11,'','minR','','maxR';...
    };
