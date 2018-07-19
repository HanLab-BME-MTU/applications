function         [idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,...
            idGroup7Selected,idGroup8Selected,idGroup9Selected] = ...
            sortIDTracks(idTracks,iGroup,outputCellArray)
    % for duplicate ids that were assigned to multiple different groups,
    % assign them to the latest group
    [uniqIdTracks, ia, ic]=unique(idTracks,'stable');
    uniqIGroup = zeros(size(uniqIdTracks));
    for jj=1:numel(ia)
        laterIdx = find(ic==jj,1,'last');
        uniqIGroup(jj)=iGroup(laterIdx);
    end
    idGroup1Selected = uniqIdTracks(uniqIGroup==1);
    idGroup2Selected = uniqIdTracks(uniqIGroup==2);
    idGroup3Selected = uniqIdTracks(uniqIGroup==3);
    idGroup4Selected = uniqIdTracks(uniqIGroup==4);
    idGroup5Selected = uniqIdTracks(uniqIGroup==5);
    idGroup6Selected = uniqIdTracks(uniqIGroup==6);
    idGroup7Selected = uniqIdTracks(uniqIGroup==7);
    idGroup8Selected = uniqIdTracks(uniqIGroup==8);
    idGroup9Selected = uniqIdTracks(uniqIGroup==9);
    if nargin>2 && outputCellArray
        idGroup1Selected = {idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                    idGroup7Selected,idGroup8Selected,idGroup9Selected};
    end
end      