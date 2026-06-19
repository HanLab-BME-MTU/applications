classdef FilopodiaTrackingProcess < DataProcessingProcess
    % Process 3. Link tip (+base) detections over time by LAP (lap.m from
    % uTrack), gap-close, then derive L(t)=geodesic base->tip, signed dL/dt
    % (protrusion +, retraction -), motility state and fluctuation frequency.
    % To reuse the full uTrack engine instead, swap the wrapper for
    % trackCloseGapsKalmanSparse on the movieInfo from process 2.

    methods (Access = public)
        function obj = FilopodiaTrackingProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                ip = inputParser;
                ip.addRequired('owner', @(x) isa(x, 'MovieData'));
                ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
                ip.addOptional('funParams', [], @isstruct);
                ip.parse(owner, varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;

                super_args{1} = owner;
                super_args{2} = FilopodiaTrackingProcess.getName;
                super_args{3} = @trackMovieFilopodia;
                if isempty(funParams)
                    funParams = FilopodiaTrackingProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end

        function varargout = loadChannelOutput(obj, iChan, varargin)
            outputList = {'adhesionTracks', 'tracksFinal'};
            ip = inputParser;
            ip.addRequired('iChan', @(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame', [], @(x) isempty(x) || all(obj.checkFrameNum(x)));
            ip.addParamValue('useCache', true, @islogical);
            ip.addParamValue('iZ', [], @(x) isempty(x) || ismember(x,1:obj.owner_.zSize_));
            ip.addParamValue('output', 'tracksFinal', @(x) all(ismember(x, outputList)));
            ip.parse(iChan, varargin{:});
            output = ip.Results.output; iFrame = ip.Results.iFrame;
            if ischar(output), output = {output}; end
            s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache);

            varargout = cell(numel(output), 1);
            for i = 1:numel(output)
                switch output{i}
                    case 'tracksFinal',    varargout{i} = s.tracksFinal;
                    case 'adhesionTracks', varargout{i} = s.adhesionTracks;
                end
                % per-frame truncation so dragtails grow as the movie plays
                if strcmp(output{i},'tracksFinal') && ~isempty(iFrame) && ~isempty(s.tracksFinal)
                    tf = s.tracksFinal;
                    sel = getTrackSEL(tf);
                    valid = (iFrame >= sel(:,1) & iFrame <= sel(:,2));
                    [tf(~valid).tracksCoordAmpCG] = deal([]);
                    nC = (iFrame - sel(valid,1) + 1) * 8;
                    vOut = tf(valid);
                    for j = 1:numel(vOut)
                        vOut(j).tracksCoordAmpCG = vOut(j).tracksCoordAmpCG(:, 1:nC(j));
                    end
                    tf(valid) = vOut;
                    varargout{i} = tf;
                end
            end
        end

        function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name = 'Tip adhesion tracks';
            output(1).var  = 'tracksFinal';
            output(1).type = 'overlay';
            output(1).formatData = @TrackingProcess.formatTracks2D;
            output(1).defaultDisplayMethod = @(x) TracksDisplay('Color', colors(x,:));
        end
    end


    methods (Access = public)
        function markSuccess(obj)
            % Mark process as successfully completed (for GUI status display).
            % success_ is SetAccess=protected so only subclass methods can set it.
            obj.success_ = true;
            obj.updated_ = true;
        end
    end

    methods (Static)
        function name = getName(), name = 'Filopodia Tracking'; end
        function h = GUI(), h = @filopodiaTrackingProcessGUI; end

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.OutputDirectory = [outputDir filesep 'FilopodiaForcePackage' filesep 'FilopodiaTracking'];
            funParams.ChannelIndex    = 1;
            funParams.DetProcessIndex = [];      % []->find FilopodiaDetectionProcess
            funParams.MaxLinkDist     = 8;       % px / frame; small: generous detection keeps tips close frame-to-frame
            funParams.LinkUseBase     = true;
            funParams.MaxGapFrames    = 3;       % frames; bridge short disappearances
            funParams.MinTrackLength  = 3;       % frames
            funParams.VelSmoothWin    = 3;       % frames
            funParams.PauseThreshVel  = [];      % length/frame; []->auto
            funParams.FreqMethod      = 'psd';   % 'psd'|'autocorr'|'cyclecount'
            funParams.DetrendL        = true;
            funParams.BatchMode       = false;
        end
    end
end
