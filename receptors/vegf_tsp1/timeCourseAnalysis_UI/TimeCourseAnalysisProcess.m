classdef TimeCourseAnalysisProcess < Process
    %MovieList Process that stores the result of timeCourseAnalysis 
    %
    %PROPERTIES
    %       summary_: Structure array storing various results for each
    %                 movie, ordered based on their time (corresponding
    %                 to caseTimeList). Contains the fields: 
    %           .diffSummary         : Diffusion analysis summary, as output
    %                                  by summarizeDiffAnRes.
    %           .diffCoefMeanPerClass: Mean diffusion coefficient per motion
    %                                  class. Order: Immobile, confined, free,
    %                                  directed, undetermined.
    %           .confRadMeanPerClass : Mean confinement radius per motion
    %                                  class. Same order as above.
    %           .ampMeanPerClass     : Mean particle amplitude per motion
    %                                  class. Rows in same order as above.
    %                                  Columns: first for "absolute" amplitude,
    %                                  second for amplitude normalized by
    %                                  monomer amplitude, as derived from modal
    %                                  analysis of particles in last 20 frames
    %                                  of each movie.
    %           .ampStatsF20         : Amplitude statistics in first 20
    %                                  frame of each movie. 
    %                                  Order: mean amplitude, first mode
    %                                  mean, first mode std, first mode
    %                                  fraction, number of modes, mean
    %                                  normalized by monomer amplitude (see
    %                                  ampMeanPerClass for details).
    %           .ampStatsL20         : Same as ampStatsF20, but for last 20
    %                                  frames of movie.

    
    properties
        summary_
    end
    %% Constructor
    methods(Access = public)
        function obj = TimeCourseAnalysisProcess(owner)
            %obj.owner_ = owner;
            %obj.name_ = getName();
            obj = obj@Process(owner, TimeCourseAnalysisProcess.getName());
            obj.summary_ = [];
        end
    end
    %% Get Set
    methods
        function setSummary(obj, summary)
            obj.summary_ = summary;
        end
    end
    %% Superclass abstracts
    methods(Static)
        function funParams = getDefaultParams()
        end
        function name = getName()
            name = 'TimeCourseAnalysisProcess';
        end
    end
end

