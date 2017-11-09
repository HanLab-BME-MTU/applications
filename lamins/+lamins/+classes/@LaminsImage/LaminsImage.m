classdef LaminsImage < matlab.mixin.SetGet
    properties
        coordinates
        steerable
        parent
        image
        adjusted
        adjustedComplement
        stretchedComplement
        skeleton
        reader
        mask
        nmsSkel
        threshSkel
        extendedSkel
        auditedSkel
        clear
        steerableFromFile = false
    end
    methods
        function obj = LaminsImage(parent,c,t,z)
            if(nargin > 0)
                if(isnumeric(parent))
                    obj.image = parent;
                    if(nargin == 4)
                        obj.coordinates = { c t z };
                    end
                elseif(isa(parent,'ProxyReader') && nargin == 1)
                    % if we are given just a ProxyReader then see 
                    obj.parent = parent;
                    proxies = parent.findProxies();
                    idx = cellfun(@(x) isa(x,'SubIndexReader'),proxies);
                    lastidx = find(idx,1,'last');
                    assert(~isempty(lastidx),'You must specify coordinates');
                    obj.coordinates = cellfun(@(x) x(1),proxies{lastidx}.subIndices);
                    obj.reader = proxies{lastidx}.reader;
                    %obj.image = proxies{lastidx}.loadImage(1,1,1);
                else
                    obj.parent = parent;
                    obj.coordinates = {c t z};
                    obj.steerableFromFile = true;
                    %obj.image;
                end
            end
        end
        function R = get.reader(obj)
            R = obj.parent.reader;
            if(~isa(R,'CellReader'))
                R = CellReader(R);
            end
        end
        function I = get.image(obj)
            if(isempty(obj.image))
                R = obj.reader;
                I = R{obj.coordinates{:}};
               obj.image = R{obj.coordinates{:}};
            end
            I = obj.image;
        end
        function A = get.adjusted(obj)
            if(isempty(obj.adjusted))
                I = obj.image;
                obj.adjusted = imadjust(obj.image,stretchlim(I,0));
            end
            A = obj.adjusted;
        end
        function C = get.adjustedComplement(obj)
            C = imcomplement(obj.adjusted);
        end
        function C = get.stretchedComplement(obj)
            C = imcomplement(imadjust(double(obj)));
        end
        function E = get.clear(obj)
            E = [];
        end
        function D = double(obj)
            D = im2double(obj.adjusted);
        end
        function steerable = get.steerable(obj)
            if(isempty(obj.steerable))
                if(~isempty(obj.parent) && ~isempty(obj.parent.params))
                    path = ['~/matlab/lamins/work/steerable_' num2str(obj.parent.params.movieNum)  '.mat'];
                else
                    obj.steerableFromFile = false;
                end
                if(obj.steerableFromFile && exist(path))
                    matobj = matfile(path);
                    obj.steerable = matobj.steerable(obj.coordinates{[1 3]});
                    obj.steerable = obj.steerable{1};
                else
                    sigma = 5;
                    if(~isempty(obj.parent) && isa(obj.parent,'lamins.classes.LaminsData'))
                        sigma = obj.parent.params.steerable.sigma;
                    end
                    disp(sigma);
                    [s.res, s.theta, s.nms] = steerableDetector(double(obj),4,sigma);
                    obj.steerable = s;
                end
            end
            steerable = obj.steerable;
        end
%         function nms = get.nms(obj)
%             S = obj.getSteerable();
%             nms = S.nms;
%         end
%         function theta = get.theta(obj)
%             S = obj.getSteerable();
%             theta = S.theta;
%         end
%         function res = getResponse(obj)
%             S = obj.getSteerable();
%             res = S.res;
%         end
        function mask = get.mask(obj)
            import lamins.functions.*;
            if(isempty(obj.mask))
                obj.mask = maskFromSteerable(obj.steerable);
            end
            mask = obj.mask;
        end
        function skel = get.nmsSkel(obj)
            if(isempty(obj.nmsSkel))
                obj.nmsSkel = obj.steerable.nms ~= 0 & obj.mask;
                obj.nmsSkel = bwmorph(obj.nmsSkel,'skel',Inf);
            end
            skel = obj.nmsSkel;
        end
        function tskel = get.threshSkel(obj)
            import lamins.functions.*;
            if(isempty(obj.threshSkel))
                theta_stats = getSegmentOrientationStats(obj.steerable.theta,obj.nmsSkel);
    %             areaThresh = thresholdRosin([theta_stats.rp.Area]);
                areaThresh = 5;
                threshed.cc = filtercc(theta_stats.cc,[theta_stats.rp.Area] > areaThresh);
                threshed.lm = labelmatrix(threshed.cc);
                obj.threshSkel = threshed.lm > 0;
            end
            tskel = obj.threshSkel;
        end
        function eskel = get.extendedSkel(obj)
            import lamins.functions.*;
            if(isempty(obj.extendedSkel))
                endpts = bwmorph(obj.threshSkel,'endpoints');
                vector = getEndPointVector(obj.threshSkel,'local');
                maxPixelExtension = 20;
                connectivity = 8;
                obj.extendedSkel = extendVectorUntilConnected(obj.threshSkel,endpts,vector,maxPixelExtension,connectivity);
            end
            eskel = obj.extendedSkel;
        end
        function askel = get.auditedSkel(obj)
            if(isempty(obj.auditedSkel))
                A = lamins.functions.auditSkelEdges(obj.extendedSkel,double(obj));
                obj.auditedSkel = A.bw;
            end
            askel = obj.auditedSkel;
        end
        function skel = get.skeleton(obj)
            import lamins.classes.*;
            if(isempty(obj.skeleton))
                obj.skeleton = Skeleton(obj);
            end
            skel = obj.skeleton;
        end
        function him = show(obj,varargin)
            him = imshow(obj.adjusted,varargin{:});
        end
        function him = imshow(obj,varargin)
            if(nargin == 1)
                varargin{1} = [];
            end
            him = imshow(obj.image,varargin{:});
        end
        function him = showMask(obj)
            if(isscalar(obj))
                him = imshow(obj.mask);
            else
                Mask = reshape({ obj.mask },size(obj));
                Mask = cell2mat(Mask);
                him = imshow(Mask);
                for i=1:size(obj,2)
                    text(512+(i-1)*1024,512,num2str(i),'color','r')
                end
            end
        end
        function out = keyboard(I)
            out = [];
            keyboard;
        end
        function him = overlayMask(obj)
            if(isscalar(obj))
                him = imshowpair(obj.adjusted,obj.mask);
            else
                Images = cell2mat( reshape({ obj.adjusted },size(obj)) );
                Mask = cell2mat( reshape({ obj.mask },size(obj)) );
                him = imshowpair(Mask,Images);
            end
        end
        function disp(obj)
            disp(['Coordinates: ']);
            obj.coordinates
            disp(['Size: ']);
%             size(obj.image)
        end
        function resetProps(obj)
            obj.image = [];
            obj.adjusted = [];
            obj.mask = [];
            obj.nmsSkel = [];
            obj.threshSkel = [];
        end
        function C = formatCell(obj,C)
            C = reshape(C,size(obj));
        end
        function clearProp(obj,prop)
            C = cell(size(obj));
            [obj.(prop)] = C{:};
        end
        function out = parget(obj,prop)
            out = cell(size(obj));
            parfor c = size(obj,1)
                for z = size(obj,2)
                    obj(c,z).(prop);
                end
            end
        end
        function [thresh,g,out] = maskThresh(obj,I)
            if(nargin < 2)
                I = obj.image;
            end
            if(all(~obj.mask(:)))
                warning('invalid mask');
                % default to 95th percentile of the image
                thresh = prctile(I(:),95);
                g = [];
                out = [];
                return;
            end
            [maskedFcn.v, maskedFcn.x] = ecdf(I(obj.mask));
            [nonmaskedFcn.v,nonmaskedFcn.x] = ecdf(I(~obj.mask));
            maskedFcn.f = @(q) interp1(double(maskedFcn.x(2:end)),maskedFcn.v(2:end),q);
            nonmaskedFcn.f = @(q) interp1(double(nonmaskedFcn.x(2:end)),nonmaskedFcn.v(2:end),q);
            g = @(q) nonmaskedFcn.f(q) - maskedFcn.f(q);
            r = linspace(double(min(I(:))),double(max(I(:))),1000);
            try
                [p,pi] = findpeaks(g(r),'sortstr','descend');
                thresh = r(pi);
                thresh = thresh(1);
            catch
                % get the 95th percentile of the unmasked region if
                % finding peaks failed
                thresh = prctile(I(~obj.mask(:)),95);
                p = NaN;
                pi = NaN;
                warning('No peaks found');
            end
            if(nargout > 2)
                out.maskedFcn = maskedFcn;
                out.nonmaskedFcn = nonmaskedFcn;
                out.maskedPct = maskedFcn.f(thresh);
                out.nonmaskedPct = nonmaskedFcn.f(thresh);
                out.p = p;
                out.pi = pi;
                out.r = r;
            end
        end
        function X = flattenIntensity(obj)
            if(isscalar(obj))
                kernel = fspecial('gaussian',50,10);
                X = double(obj)./imfilter(double(obj),kernel);
            else
                X = arrayfun(@flattenIntensity,obj,'UniformOutput',false);
            end
        end
        function s = saveSteerable(obj)
            path = ['~/matlab/lamins/work/steerable_' num2str(obj(1).parent.params.movieNum)  '.mat'];
            assert(~exist(path),[path ' exists. Not saving steerable']);            
            S.steerable = {obj.steerable};
            save(path,'-struct','S','-v7.3');
            s = true
        end
        function loadSteerable(obj)
            path = ['~/matlab/lamins/work/steerable_' num2str(obj(1).parent.params.movieNum)  '.mat'];
            if(exist(path))
                S = load(['~/matlab/lamins/work/steerable_' num2str(obj(1).parent.params.movieNum)  '.mat']);
                [obj.steerable] = S.steerable{:};
            end
        end
        function c = getMaskCircularity(obj)
            % calculate nucleus shape factor from mask
            % a perfect circle has a circularity of 1
            % a straight line has a circularity of 0
            rp = regionprops(obj.mask,'Area','Perimeter');
            assert(isscalar(rp));
            c = 4*pi*(rp.Area)/rp.Perimeter.^2;
        end
    end
end
