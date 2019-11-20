function [] = exportTrackAmpsFromTracksNA(MD,iChan)
%function [] = exportTrackAmpsFromTracksNA(MD) reads G1 and G2 tracks from
%focal adhesion package and export only amp and amp2 (for now)
%   input
%       MD:     movieData, where focal adhesion package has been run.
%       iChan:  the channel index of one analyzed
%   output
%       It will create a folder and export tracks variable with amp and
%       amp2 as members
% Sangyoon Han, November 2019

% Load FA package
faPack = MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load idsClasses
iInitialRise = 11; iClassProc = 8;
initRiseProc = faPack.getProcess(iInitialRise);
classProc = faPack.getProcess(iClassProc);
% chanOutput = initRiseProc.checkChannelOutput;
chanOutput = classProc.checkChannelOutput;
if sum(chanOutput)>0
    iChan=find(chanOutput);
end
idsClassfied = classProc.loadChannelOutput(iChan,'output','iClassesAll');
% idsClassfied = initRiseProc.loadChannelOutput(iChan,'output','idClass');

% % test
% trackslong = initRiseProc.loadChannelOutput(iChan,'output','tracksNA');

% Load tracks selectively
tracksG1long = initRiseProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',find(idsClassfied.idGroup1));
tracksG2long = initRiseProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',find(idsClassfied.idGroup2));
% Maybee some inspection here

% Export to simpler variables
myFieldNames={'amp','amp2','ampTotal','ampTotal2'};
nElementsG1 = numel(tracksG1long);
nElementsG2 = numel(tracksG2long);
tracksG1(nElementsG1,1) = struct(myFieldNames{1},[],myFieldNames{2},[],myFieldNames{3},[],myFieldNames{4},[]);
tracksG2(nElementsG2,1) = struct(myFieldNames{1},[],myFieldNames{2},[],myFieldNames{3},[],myFieldNames{4},[]);

% try
    for curFN=myFieldNames
        for ii=nElementsG1:-1:1
            try
                tracksG1(ii,1).(curFN{1}) = tracksG1long(ii,1).(curFN{1});
            catch
                disp(['Skipping ' num2str(ii) 'th tracks'])
            end
        end
        for ii=nElementsG2:-1:1
            try
                tracksG2(ii,1).(curFN{1}) = tracksG2long(ii,1).(curFN{1});
            catch
                disp(['Skipping ' num2str(ii) 'th tracks'])
            end
        end
    end
% catch
%     myFieldNames2={'ampTotal','ampTotal2'};
%     for curFN=myFieldNames
%         for ii=nElementsG1:-1:1
%             tracksG1(ii,1).(curFN{1}) = tracksG1long(ii,1).(curFN{1});
%         end
%         for ii=nElementsG2:-1:1
%             tracksG2(ii,1).(curFN{1}) = tracksG2long(ii,1).(curFN{1});
%         end
%     end
% end

% Saving
faDir = [faPack.outputDirectory_ filesep 'TracksSimple'];
g1Dir = [faDir filesep 'G1'];
g2Dir = [faDir filesep 'G2'];
if ~exist(g1Dir,'dir')
    mkdir(g1Dir); mkdir(g2Dir);
end

save([g1Dir filesep 'tracksG1.mat'], 'tracksG1')
save([g2Dir filesep 'tracksG2.mat'], 'tracksG2')

disp('done!')
end

