classdef OrientationScaleMatrix < double
    %OrientationScaleMatrix Matrix that represents orientation x scale
    %
    % Rows represent orientation and are equispaced samples of a fourier
    % series
    %
    % Cols represent scale on an interval and are sampled on a
    % Chebyshev-Lobatto grid
    %
    % () indexing will retrieve the underlying doubles or if appropriate
    % another OrientationScaleMatrix
    % An OrientationScaleMatrix is returned if a two element vector is
    % indexed in the scale dimension
    %
    % {} indexing attempts to convert the rows and columns into chebfuns
    %
    % This is a subclass of double, so the original values can be retrieved
    % by casting back to a double. Except for the above
    %
    % transposition is ignored except for
    % 1) Use in imagesc to rotate the plot
    % 2) Use in max to determine default dimension
    
    properties
        scaleRange;
        orientationRange;
        rotated;
    end
    
    methods
        function obj = OrientationScaleMatrix(matrix,scaleRange,orientationRange,rotated)
            obj = obj@double(matrix);
            if(nargin > 1)
                obj.scaleRange = scaleRange;
            else
                obj.scaleRange = [0 1];
            end
            if(nargin > 2)
                obj.orientationRange = orientationRange;
            else
                obj.orientationRange = [0 360];
            end
            if(nargin > 3)
                obj.rotated = rotated;
            else
                obj.rotated = false;
            end
        end
        function varargout = imagesc(obj,varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.addParameter('rotate',obj.rotated,@islogical);
            ip.parse(varargin{:});
            E = obj.sample(360,101);
            scaleSeq = linspace(obj.scaleRange(1),obj.scaleRange(2),size(E,2));
            orientationSeq = linspace(obj.orientationRange(1),obj.orientationRange(2) - diff(obj.orientationRange)/size(E,1),size(E,1));
            if(~ip.Results.rotate)
                [varargout{1:nargout}] = imagesc(scaleSeq,orientationSeq,E,ip.Unmatched);
                xlabel('Scale');
                ylabel('Orientation (Degrees)');
                grid on;
                set(gca,'YTick',[0 obj.orientationRange(2)/6:obj.orientationRange(2)/6:obj.orientationRange(2)- diff(obj.orientationRange)/360]);
                set(gca,'YMinorTick','on');
                set(gca,'YMinorGrid','on');
            else
                [varargout{1:nargout}] = imagesc(orientationSeq,scaleSeq,E.',ip.Unmatched);
                ylabel('Scale');
                xlabel('Orientation (Degrees)');
                grid on;
                set(gca,'XTick',[0 obj.orientationRange(2)/6:obj.orientationRange(2)/6:obj.orientationRange(2)- diff(obj.orientationRange)/360]);
                set(gca,'XMinorTick','on');
                set(gca,'XMinorGrid','on');
            end
        end
        function varargout = contour(obj,varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.addParameter('rotate',obj.rotated,@islogical);
            ip.parse(varargin{:});
            E = obj.sample(360,101);
            scaleSeq = linspace(obj.scaleRange(1),obj.scaleRange(2),size(E,2));
            orientationSeq = linspace(obj.orientationRange(1),obj.orientationRange(2) - diff(obj.orientationRange)/size(E,1),size(E,1));
            if(~ip.Results.rotate)
                [varargout{1:nargout}] = contour(scaleSeq,orientationSeq,E);
                xlabel('Scale');
                ylabel('Orientation (Degrees)');
                grid on;
                set(gca,'YTick',[0 obj.orientationRange(2)/6:obj.orientationRange(2)/6:obj.orientationRange(2)- diff(obj.orientationRange)/360]);
                set(gca,'YMinorTick','on');
                set(gca,'YMinorGrid','on');
            else
                [varargout{1:nargout}] = contour(orientationSeq,scaleSeq,E.',ip.Unmatched);
                ylabel('Scale');
                xlabel('Orientation (Degrees)');
                grid on;
                set(gca,'XTick',[0 obj.orientationRange(2)/6:obj.orientationRange(2)/6:obj.orientationRange(2)- diff(obj.orientationRange)/360]);
                set(gca,'XMinorTick','on');
                set(gca,'XMinorGrid','on');
            end
        end
        function B = transpose(A)
            B = A;
            B.rotated = true;
        end
        function B = ctranspose(A)
            B = transpose(A);
        end
        function E = expand(obj,nOrientations,nScales)
            E = double(obj);
            
            if(nargin > 2 && nScales > size(obj,2))
                E = [E E(:,end-1:-1:2)];
                E = interpft1([-pi pi],E.',linspace(pi,0,nScales).','horner').';
                E = E(:,(end-nScales+1):end);
            end
            if(~isempty(nOrientations) && nOrientations > size(obj,1))
                E = interpft(E,nOrientations,1);
            end
        end
        function E = sample(obj,nOrientations,nScales)
            %sample at equally spaced points in actual space
            E = double(obj);
            
            if(nargin > 2 && nScales > size(obj,2))
                E = [E E(:,end-1:-1:2)];
                E = interpft1([-pi pi],E.',acos(linspace(-1,1,nScales)).','horner').';
                E = E(:,(end-nScales+1):end);
            end
            if(~isempty(nOrientations) && nOrientations > size(obj,1))
                E = interpft(E,nOrientations,1);
            end
        end
        function B = subsref(A,S)
            outScaleRange = A.scaleRange;
            switch(S(1).type)
                case '()'
                    assert(length(S(1).subs) == 2);
                    recast = true;
                    % Do nothing if : is given
                    if(S(1).subs{1} == ':')
                        B = double(A);
                    else
                        B = interpft1(A.orientationRange,double(A),S(1).subs{1}(:));
                        recast = false;
                    end
                    if(S(1).subs{2} == ':')
                    else
                        
                        if(length(S(1).subs{2}) == 2)
                            outScaleRange = S(1).subs{2};
                            S(1).subs{2} = chebpts(size(A,2),S(1).subs{2}).';
                        else
                            recast = false;
                        end
                        chebAbscissa = acos((S(1).subs{2} - A.scaleRange(1))/(A.scaleRange(2)-A.scaleRange(1))*2-1);
                        B = [B B(:,end-1:-1:2)];
                        B = interpft1([-pi pi],B.',chebAbscissa(:),'horner').';
                    end
                    if(recast && ~isa(B,'OrientationScaleMatrix'))
                        B = OrientationScaleMatrix(B,outScaleRange,A.orientationRange);
                    end
                case '{}'
                    assert(length(S(1).subs) == 2);
%                     assert(isscalar(S(1).subs{1}) && isscalar(S(1).subs{2}));
                    S(1).type = '()';
                    B = subsref(A,S(1));
                    if(strncmp(class(B),'double',6))
                        if(isrow(B) || ...
                                isscalar(S(1).subs{2}) && S(1).subs{2} == ':')
                            B = chebfun(B.',A.scaleRange);
                        elseif(iscolumn(B) || ...
                                isscalar(S(1).subs{1}) && S(1).subs{1} == ':')
                            B = chebfun(B,A.orientationRange,'trig');
                        elseif(length(S(1).subs{2}) == 2)
                            B = chebfun(B.',S(1).subs{2});
                        end
                    elseif(isa(B,'OrientationScaleMatrix'))
                        B = double(B);
                        B = chebfun2([B B(:,end-1:-1:2)],[-1 1 A.orientationRange],'trig');
                    end
                otherwise
                    B = subsref@double(A,S);
                    if(strncmp(class(B),'double',6) && ~isvector(B))
                        try
                            B = OrientationScaleMatrix(B,A.scaleRange,A.orientationRange);
                        catch err
                        end
                    end
                    return;
            end
            if(length(S) > 1)
                B = subsref(B,S(2:end));
            end
        end
        function varargout = max(A,B,dim,varargin)
            if(nargin < 2)
                B = [];
            end
            if(nargin < 3)
                dim = 1 + A.rotated;
            end
            assert(dim == 1 || dim == 2);
            
            if(~isvector(B) && ~isempty(B))
                [varargout{1:nargout}] = max@double(A,B,dim,varargin{:});
                return;
            end
            
            if(dim == 1)
                [varargout{1:nargout}] = max(chebfun(double(A),A.orientationRange,'trig'),varargin{:});
            elseif(dim == 2)
                [varargout{1:nargout}] = max(chebfun(double(A).',A.scaleRange),varargin{:});
                if(A.rotated)
                    varargout{1} = varargout{1}.';
                    varargout{2} = varargout{2}.';
                end
            end
            
            if(~isempty(B))
                varargout{1} = max(varargout{1},B);
                varargout = varargout(1);
            end
        end
        function varargout = min(A,B,dim,varargin)
            if(nargin < 2)
                B = [];
            end
            if(nargin < 3)
                dim = 1 + A.rotated;
            end
            assert(dim == 1 || dim == 2);
            
            if(~isvector(B) && ~isempty(B))
                [varargout{1:nargout}] = min@double(A,B,dim,varargin{:});
                return;
            end
            
            if(dim == 1)
                [varargout{1:nargout}] = min(chebfun(double(A),A.orientationRange,'trig'),varargin{:});
            elseif(dim == 2)
                [varargout{1:nargout}] = min(chebfun(double(A).',A.scaleRange),varargin{:});
                if(A.rotated)
                    varargout{1} = varargout{1}.';
                    varargout{2} = varargout{2}.';
                end
            end
            
            if(~isempty(B))
                varargout{1} = min(varargout{1},B);
                varargout = varargout(1);
            end
        end
    end   
end