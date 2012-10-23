function corrMoviesEBtoKinDynamics(movieDataOrList,varargin)
%CORRMOVIESEBTOKINDYNAMICS correlates EB signal at kinetochores to kinetochore dynamics
%
%Khuloud Jaqaman, 10/2012

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieDataOrList', @(x) (isa(x,'MovieData')||isa(x,'MovieList')));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrList,varargin{:});
paramsIn=ip.Results.paramsIn;

%determine whether what is input is of type MovieData or MovieList
if isa(movieDataOrList,'MovieList')
    inputList = 1;
    numCell = length(movieDataOrList.movieDataFile_);
else
    inputList = 0;
    numCell = 1;
end

%Get the indices of any previous EB-kin correlation processes
%If the process doesn't exist, create it
iProc = movieDataOrList.getProcessIndex('CorrEBtoKinDynamicsProcess',1,0);
if isempty(iProc)
    iProc = numel(movieDataOrList.processes_)+1;
    movieDataOrList.addProcess(CorrEBtoKinDynamicsProcess(movieDataOrList));
end
corrProc = movieDataOrList.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(corrProc,paramsIn);

%% --------------- Initialization ---------------%%

% Check EB detection process and spindle pole process
% load their output
% also set up the input file paths for this process
sisterListEB = cell(numCell,1);
poleInfo = cell(numCell,1);
ebChan = NaN(numCell,1);
inFilePaths = cell(numCell,2);
if inputList
    
    for iCell = 1 : numCell
        
        movieData = MovieData.load(movieDataOrList.movieDataFile_{iCell});
        
        %EB detection
        iEBDetProc = movieData.getProcessIndex('KEBDetectionProcess',1,1);
        
        assert(~isempty(iEBDetProc),['EB detection has not been run for Cell ' num2str(iCell) '! '...
            'Please run EB detection prior to correlating it with kinetochore dynamics!'])
        ebDetProc = movieData.processes_{iEBDetProc};
        
        ebChan(iCell) = find(ebDetProc.checkChannelOutput);
        
        assert(all(ebDetProc.checkChannelOutput(ebChan(iCell))),...
            ['Missing EB detection output for Cell ' num2str(iCell) '! Please run EB detection before ' ...
            'correlating it with kinetochore dynamics!']);
        
        sisterListEB{iCell} = ebDetProc.loadChannelOutput(ebChan(iCell));
        inFilePaths{iCell,1} = ebDetProc.outFilePaths_{1,ebChan(iCell)};

        %spindle pole
        iPoleProc = movieData.getProcessIndex('SpindlePolesEBProcess',1,1);
        
        assert(~isempty(iPoleProc),['Spindle poles have not been detected for Cell ' num2str(iCell) '! '...
            'Please run spindle pole detection prior to correlating EB signal with kinetochore dynamics!'])
        poleProc = movieData.processes_{iPoleProc};
        
        assert(all(poleProc.checkChannelOutput(ebChan(iCell))),...
            ['Missing spindle pole output for Cell ' num2str(iCell) '! Please detect spindle poles before ' ...
            'correlating EB signal with kinetochore dynamics!']);
        
        poleInfo{iCell} = poleProc.loadChannelOutput(ebChan(iCell));
        inFilePaths{iCell,2} = poleProc.outFilePaths_{1,ebChan(iCell)};
        
    end
    
else
    
    %EB detection
    iEBDetProc = movieDataOrList.getProcessIndex('KEBDetectionProcess',1,1);
    
    assert(~isempty(iEBDetProc),['EB detection has not been run! '...
        'Please run EB detection prior to correlating it with kinetochore dynamics!'])
    ebDetProc = movieDataOrList.processes_{iEBDetProc};
    
    ebChan = find(ebDetProc.checkChannelOutput);
    
    assert(all(ebDetProc.checkChannelOutput(ebChan)),...
        ['Missing EB detection output! Please run EB detection before ' ...
        'correlating it with kinetochore dynamics!']);
    
    sisterListEB{1} = ebDetProc.loadChannelOutput(ebChan);
    inFilePaths{1,1} = ebDetProc.outFilePaths_{1,ebChan};
    
    %spindle pole
    iPoleProc = movieDataOrList.getProcessIndex('SpindlePolesEBProcess',1,1);
    
    assert(~isempty(iPoleProc),['Spindle poles have not been detected! '...
        'Please run spindle pole detection prior to correlating EB signal with kinetochore dynamics!'])
    poleProc = movieDataOrList.processes_{iPoleProc};
    
    assert(all(poleProc.checkChannelOutput(ebChan)),...
        ['Missing spindle pole output! Please detect spindle poles before ' ...
        'correlating EB signal with kinetochore dynamics!']);
    
    poleInfo{1} = poleProc.loadChannelOutput(ebChan);
    inFilePaths{1,2} = poleProc.outFilePaths_{1,ebChan};
    
end
corrProc.setInFilePaths(inFilePaths);

% Set up the output file
outFilePaths{1} = [p.OutputDirectory filesep 'channels_12.mat'];
mkClrDir(p.OutputDirectory);
corrProc.setOutFilePaths(outFilePaths);

%% --------------- EB - kinetochore dynamics correlation ---------------%%%

disp('Correlating EB & kinetochore dynamics...')

%call function
measurementsEB = corrEBtoKinDynamics(sisterListEB,poleInfo,'minDisp',p.minDisp); %#ok<NASGU>

%save output
save(outFilePaths{1},'measurementsEB')

disp('Done')
