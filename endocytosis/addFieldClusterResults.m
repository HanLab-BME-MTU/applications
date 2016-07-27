function [experiment] = addFieldClusterResults(experiment);

%save old directory
oldDir = cd;

%FOR EACH MOVIE
for iexp = 1:length(experiment)

    waitHandle = waitbar(iexp/length(experiment),['running movie ' num2str(iexp) ' out of ' num2str(length(experiment))]);
    
    %FIND CLUSTER RESULSTS
    %go to movie folder
    cd(experiment(iexp).source)

    if exist('ClusterData','dir') == 7
        %move to clustering results folder
        cd('ClusterData')
        %find all files under folder
        files = dir('*.mat');

        %find the ones that contain clusterResults in the name
        findClusterResults = [];
        for ifile = 1:length(files)
            findClusterResults(ifile) = ~isempty(strfind(files(ifile).name,'lusterResults'));
        end
        clusterResultsIndex = find(findClusterResults == 1);
        %if more than one clusterResult ask the user to pick wihc one he wants
        %to use
        if length(clusterResultsIndex) > 1
            %list all movies and ask user to select movies to use in analysis
            [selection, selectionList] = listSelectGUI({files(clusterResultsIndex).name},1,'move',[]);
            filesName = selectionList{1};
            load(filesName)
            experiment(iexp).clusterResults = clusterResults;
        elseif length(clusterResultsIndex) == 1
            files = files(clusterResultsIndex);
            load(files.name)
            experiment(iexp).clusterResults = clusterResults;
        else
            experiment(iexp).clusterResults = [];
        end

        %in the future it might be nice if the user gets to pick which ones to
        %add; although, it might be easier to just add the latest ones and
        %change individually the ones that are not the latest
    else
        experiment(iexp).clusterResults = [];
    end %of if clusterResults found

    close(waitHandle);
    
end %of each movie



%return to old directory
cd(oldDir)

end %of function