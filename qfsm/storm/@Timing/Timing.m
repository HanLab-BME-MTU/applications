classdef Timing < handle
    
    % ---------------------
    % This class provides tools to analyze timings of an algorithm. Keep in
    % mind that the measured times might be inaccurate. A Timing object can 
    % control several timers which can be identified by a tag. Furthermore
    % data can be appended to timers.
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties (GetAccess = 'public',SetAccess = 'private')
        ticID; % The ID of the initial tic
        times; % Columns: tStart, tStop, tag, tDelta, data
    end
    
    methods
        function obj = Timing()
            obj.ticID = tic;
            obj.times = cell(0,5);
        end
        
        start(obj,tag);
        tDelta = stop(obj,tag,data);
        display(obj,tag);
        tDelta = get(obj,tag);
        tDeltaStr = getFormatted(obj,tag);
        data = getData(obj,tag);
        running(obj);
        save(obj,fullPath);
    end
    
    methods(Static = true)
        obj = load(fullPath);
    end
    
end