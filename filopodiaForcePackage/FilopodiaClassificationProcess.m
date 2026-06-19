classdef FilopodiaClassificationProcess < DataProcessingProcess
    % Process 4. Assemble filopodia from tracked tip adhesions by tracing each
    % shaft from the tip to the body along the steerable ridge orientation.
    % The trace endpoint on the body is the base; shaft arc length is the
    % filopodium length L(t). Tracks whose shaft does not reach the body often
    % enough (background blobs) are rejected. No ridge-mask skeletons are used,
    % so per-frame tips are stable (they are tracked intensity blobs).

    methods (Access = public)
        function obj = FilopodiaClassificationProcess(owner, varargin)
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
                super_args{2} = FilopodiaClassificationProcess.getName;
                super_args{3} = @classifyMovieFilopodia;
                if isempty(funParams)
                    funParams = FilopodiaClassificationProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end

        function varargout = loadChannelOutput(obj, iChan, varargin)
            outputList = {'filopodia', 'tipTracks', 'tipPos', 'basePos', 'shaftPos', 'shaftLines', 'roleByTrack'};
            ip = inputParser;
            ip.addRequired('iChan', @(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame', [], @(x) isempty(x) || all(obj.checkFrameNum(x)));
            ip.addParamValue('useCache', true, @islogical);
            ip.addParamValue('iZ', [], @(x) isempty(x) || ismember(x,1:obj.owner_.zSize_));
            ip.addParamValue('output', 'tipPos', @(x) all(ismember(x, outputList)));
            ip.parse(iChan, varargin{:});
            output = ip.Results.output; iFrame = ip.Results.iFrame;
            if ischar(output), output = {output}; end
            if isempty(iFrame), iFrame = 1; end
            s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache);

            varargout = cell(numel(output), 1);
            for i = 1:numel(output)
                switch output{i}
                    case 'filopodia',   varargout{i} = s.filopodia;
                    case 'tipTracks'
                        tf = s.tipTracks;
                        if ~isempty(tf) && ~isempty(ip.Results.iFrame)
                            sel = getTrackSEL(tf);
                            valid = (iFrame >= sel(:,1) & iFrame <= sel(:,2));
                            [tf(~valid).tracksCoordAmpCG] = deal([]);
                            nC = (iFrame - sel(valid,1) + 1) * 8;
                            vOut = tf(valid);
                            for q = 1:numel(vOut)
                                vOut(q).tracksCoordAmpCG = vOut(q).tracksCoordAmpCG(:, 1:nC(q));
                            end
                            tf(valid) = vOut;
                        end
                        varargout{i} = tf;
                    case 'roleByTrack', varargout{i} = s.roleByTrack;
                    case 'tipPos',      varargout{i} = s.posByFrame.tip{iFrame};
                    case 'basePos',     varargout{i} = s.posByFrame.base{iFrame};
                    case 'shaftPos',    varargout{i} = s.posByFrame.shaft{iFrame};
                    case 'shaftLines',  varargout{i} = s.shaftByFrame{iFrame};
                end
            end
        end

        function output = getDrawableOutput(obj)
            output(1).name = 'Filopodium shafts';
            output(1).var  = 'shaftLines';
            output(1).type = 'overlay';
            output(1).formatData = [];
            output(1).defaultDisplayMethod = @(x) LineDisplay('Color',[1 1 0],'LineWidth',1);
            output(2).name = 'Tip adhesions';
            output(2).var  = 'tipPos';
            output(2).type = 'overlay';
            output(2).formatData = [];
            output(2).defaultDisplayMethod = @(x) LineDisplay('Marker','o', ...
                'LineStyle','none','Color',[1 0 0],'MarkerSize',9);
            output(3).name = 'Shaft adhesions';
            output(3).var  = 'shaftPos';
            output(3).type = 'overlay';
            output(3).formatData = [];
            output(3).defaultDisplayMethod = @(x) LineDisplay('Marker','.', ...
                'LineStyle','none','Color',[1 1 0],'MarkerSize',8);
            output(4).name = 'Base adhesions';
            output(4).var  = 'basePos';
            output(4).type = 'overlay';
            output(4).formatData = [];
            output(4).defaultDisplayMethod = @(x) LineDisplay('Marker','o', ...
                'LineStyle','none','Color',[0 1 1],'MarkerSize',7);
            output(5).name = 'Tip tracks';
            output(5).var  = 'tipTracks';
            output(5).type = 'overlay';
            output(5).formatData = @TrackingProcess.formatTracks2D;
            output(5).defaultDisplayMethod = @(x) TracksDisplay('Color',[1 0.3 0]);
        end
    end

    methods (Static)
        function name = getName(), name = 'Filopodia Classification'; end
        function h = GUI(), h = @filopodiaClassificationProcessGUI; end

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.ChannelIndex    = 1;
            funParams.TrkProcessIndex = [];     % []->find FilopodiaTrackingProcess
            funParams.SegProcessIndex = [];     % []->find FilopodiaSegmentationProcess
            funParams.DetProcessIndex = [];     % []->find FilopodiaDetectionProcess
            % --- tip eligibility (only WELL-TRACKED adhesions drive geometry) ---
            funParams.MinTipLifetime  = 5;      % frames; tip track must persist this long
            funParams.MinLinearFrac   = 0.85;   % trajectory variance fraction on principal axis (linear)
            funParams.MinTrajSpread   = 3;      % px; below this a track is "stationary" (linearity n/a)
            funParams.MinTipDist      = 6;      % px; tip must reach at least this far from body
            % --- shaft assignment & joint assembly ---
            funParams.ShaftBand       = 4;      % px; adhesion within this of tip->base line = shaft adhesion
            funParams.WShaft          = 0.0;    % 0 = steerable only; raise to reward passing shaft adhesions (use as tie-breaker)
            funParams.WLen            = 0.25;   % small length penalty (smaller -> deeper base)
            funParams.WPrior          = 0.0;    % 0 = no neighbor prior; raise for smoother angular field across filopodia
            funParams.WOverlap        = 0.8;    % penalty if shaft crosses an already-fixed shaft
            funParams.WBaseSep        = 0.7;    % penalty if base is near an already-fixed base
            funParams.MinBaseSep      = 8;      % px; minimum spacing between distinct filopodium bases
            funParams.NeighRadius     = 60;     % px; neighborhood for direction prior
            % --- straight shaft direction (see assembleFilopodiaFrame) ---
            funParams.MaxShaftLen     = 160;    % px; max shaft length (tip reach)
            funParams.SweepRange      = 35;     % deg; direction sweep around ridge/body-ward
            funParams.SweepStep       = 3;      % deg; sweep step
            funParams.BodyMaxAngle    = 75;     % deg; shaft must be within this of body-ward (wider for tangential filopodia)
            funParams.AlignBand       = 3.5;    % px; perpendicular band for collinear tip support
            funParams.AlignWeight     = 0.12;   % weight of collinear support in direction score
            funParams.LenPenalty      = 0.6;    % penalty on shaft length (prefer radial)
            funParams.MinReachFrac    = 0.5;    % accept tip track if it acts as tip in >= this fraction of its frames
            funParams.VelSmoothWin    = 3;      % frames; smoothing for dL/dt
            funParams.OutputDirectory = [outputDir filesep 'FilopodiaForcePackage' filesep 'FilopodiaClassification'];
        end
    end
end
