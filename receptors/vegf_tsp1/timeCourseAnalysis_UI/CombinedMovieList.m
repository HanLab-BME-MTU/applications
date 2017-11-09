classdef CombinedMovieList < MovieObject
    %List of MovieLists used for timeCourseAnalysis. All ML in CML should be under same conditions
    %
    %PROPERTIES
    %   fileName_       : The name under which the combined movie list is saved
    %   filePath_       : The path where the combined movie list is saved
    %   analysisPara_   : Parameters to be used when doing timeCourseAnalysis
    %       .alignEvent     : string describing event to which ML times
    %                         should be aligned to using TimePoints process
    %   name_           : (String) Name (or experimental condition) of the dataset. Used for plotting, etc.
    %
    %Tae H Kim, July 2015
    
    properties
        fileName_           % The name under which the combined movie list is saved
        filePath_           % The path where the combined movie list is saved
        analysisPara_   %: Parameters to be used when doing timeCourseAnalysis
        name_               % Name (or experimental condition) of the dataset. Used for plotting, etc.
        %alignIndx_          % array of index numbers pointing to MDs to be used as the point of reference for aligning ML times. ie). if alginIndx_(2) = 3, the second ML is aligned using third MD in that movieList. If 0, aligned using default relativeZero
        %relTime_            % array of relative time of the point of reference MD. ie). if relTime_(2) = 4.5, the second ML's point of reference MD time is at 4.5 min realtive to other MLs.
    end
    
    properties(SetAccess = protected)
        movieListDirectory_ % Cell array of ML directory
    end
    
    properties(Transient = true)
        movieLists_         % Array of movieLists
    end
    %% Constructor
    methods
        %INPUT
        %   movieLists : array of MovieList or cell array of strings of
        %                full filepath
        function obj = CombinedMovieList(movieLists, outputDirectory, varargin)
            if isa(movieLists, 'MovieList')
                obj.movieListDirectory_ = arrayfun(@(x) x.getFullPath(), movieLists, 'UniformOutput', false);
            elseif iscellstr(movieLists)
                obj.movieListDirectory_ = movieLists;
            else
                error('invalid movieLists')
            end
            obj.createTime_ = clock;
            %nML = numel(obj.movieListDirectory_);
            %input parser to assign alignIndx, relTime, and relativeZero
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched = true;
            ip.addParameter('name', '', @(x) ischar(x) || iscellstr(x));
            ip.parse(varargin{:});
            obj.outputDirectory_ = outputDirectory;
            obj.name_ = ip.Results.name;
            %set default parameter values
            obj.analysisPara_ = CombinedMovieList.getDefaultPara();
        end
    end
    methods(Static)
        %Default values for analysisPara_
        function analysisPara = getDefaultPara()
            analysisPara.alignEvent = 'start'; %or one could also use 'VEGF_added'
        end
    end
    %% Set Get
    methods
    end
    %% Add and Remove ML
    methods
        %add
        function addML(obj, movieList, varargin)
            nML = numel(obj.movieListDirectory_) + 1;
            if isa(movieList, 'MovieList')
                obj.movieListDirectory_{nML} = movieList.getFullPath();
            elseif ischar(movieList)
                obj.movieListDirectory_{nML} = movieList;
            else
                error('invalid movieList')
            end
            %obj.alignIndx_(nML) = 0;
            %obj.relTime_(nML) = 0;
            %obj.relativeTimeZero_{nML} = clock();
        end
        %delete
        function deleteML(obj, index)
            mask = true(1,numel(obj.movieListDirectory_));
            mask(index) = false;
            obj.movieListDirectory_ = obj.movieListDirectory_(mask);
            %obj.alignIndx_ = obj.alignIndx_(mask);
            %obj.relTime_ = obj.relTime_(mask);
            %obj.relativeTimeZero_ = obj.relativeTimeZero_(mask);
        end
    end
    %% Loading, Saving, and SanityCheck
    methods
        %Save
        function save(CML, varargin)
            % Check path validity for movie list
            fullPath = CML.getFullPath();
            assert(~isempty(fullPath), 'Invalid path');
            % Save
            save(fullPath, 'CML');
        end
        %Sanity check
        function sanityCheck(obj, varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched = true;
            ip.addParameter('suppressPrinting', true, @islogical);
            ip.parse(varargin{:});
            suppressPrinting = ip.Results.suppressPrinting;
            %gets rid of suppressPrinting from the varargin. So it can be
            %passed into super.sanityCheck
            inputIndx = find(cellfun(@(x) strcmpi(x, 'suppressPrinting'), varargin));
            mask = true(1,numel(varargin));
            if ~isempty(inputIndx)
                mask([inputIndx, inputIndx+1]) = [false false];
                varargin = varargin(mask);
            end
            % Call the superclass sanityCheck
            sanityCheck@MovieObject(obj, varargin{:});
            nML = numel(obj.movieListDirectory_);
            progressTextMultiple('Loading Movie Lists', nML)
            if suppressPrinting
                for iML = 1:nML
                    evalc('ML(iML) = MovieList.load(obj.movieListDirectory_{iML});');
                    progressTextMultiple();
                end
            else
                for iML = 1:nML
                    ML(iML) = MovieList.load(obj.movieListDirectory_{iML});
                    progressTextMultiple();
                end
            end
            obj.movieLists_ = ML;
        end
    end
    methods(Static)
        %Load
        function obj = load(fullPath)
            assert(strcmpi(fullPath(end-3:end), '.mat'), 'Input must be a MAT file');
            obj = MovieObject.loadMatFile('CombinedMovieList', fullPath);
            obj.sanityCheck();
        end
    end
    %% Superclass Abstract
    methods(Static)
        %MovieObject Abstract
        function propName = getPathProperty()
            propName = 'filePath_';
        end
        %MovieObject Abstract
        function propName = getFilenameProperty()
            propName = 'fileName_';
        end
        %????????????????????????
        function status=checkValue(property,value)
           % Return true/false if the value for a given property is valid
            
           % Parse input
           ip = inputParser;
           ip.addRequired('property',@(x) ischar(x) || iscell(x));
           ip.parse(property);

           
           % Return result of validation
           status = true;
        end
    end
end

