classdef OrientationScaleSpaceWaveletFilter < OrientationSpaceFilter
    %OrientationScaleSpaceFilter implements a polar separable filter in
    %frequency space using an angular Gaussian for orientation and a
    %raised cos log filter (a la Simoncelli) for scale
    
    properties
    end
    
    methods
        function obj = OrientationScaleSpaceFilter(f_c,b_f,K,normEnergy)
            if(nargin < 3)
                K = b_f;
                b_f = NaN;
            end
            if(nargin < 4 && ischar(K))
                normEnergy = K;
                K = b_f;
                b_f = NaN;
            end
            if(nargin < 4 || isempty(normEnergy))
                normEnergy = 'none';
            else
                if(iscell(normEnergy))
                    normEnergy = normEnergy{1};
                end
            end
            obj = obj@OrientationSpaceFilter(f_c, b_f, K, normEnergy);
        end
        function setupFilter(obj,siz)
            if(isscalar(siz))
                siz = siz([1 1]);
            end
            coords = orientationSpace.getFrequencySpaceCoordinates(siz);
            notSetup = ~cellfun(@(x) isequal(siz,x),{obj.size});
            notSetup = notSetup | cellfun('isempty',{obj.F});
            obj = obj(notSetup);
            if(isempty(obj))
                return;
            end
            [obj.size] = deal(siz);
            if( all(obj(1).K == [obj.K]) )
                % angular component is all the same
                
                A = orientationSpace.angularKernel(obj(1).K,obj(1).angles,coords);
                R = radialWaveletKernel([obj.f_c], coords);
                for o=1:numel(obj)
                    obj(o).F = bsxfun(@times,A, R(:,:,o));
                end
            else
                for o=1:numel(obj)
                    obj(o).setupFilter;
%                     obj(o).F = orientationSpace.kernel(obj(o).f_c, obj(o).b_f, obj(o).K, obj(o).angles, coords);
                end
            end
%             for o=1:numel(obj)
%                 if(isempty(obj(o).normEnergy))
%                     break;
%                 end
%                 switch(obj(o).normEnergy)
%                     case 'energy'
%                         % E is complex
%                         E = shiftdim(obj(o).getEnergy(),-1);
%                         F = obj(o).F;
%                         obj(o).F = bsxfun(@rdivide,real(F),real(E)) +1j*bsxfun(@rdivide,imag(F),imag(E));
%                     case 'peak'
%                         F = obj(o).F;
%                         sumF = sum(F(:))./numel(F);
%                         obj(o).F = real(F)./real(sumF) + 1j*imag(F)./imag(sumF);
%                     case 'scale'
%                         obj(o).F = obj(o).F ./ obj(o).f_c ./ sqrt(siz(1)*siz(2));
%                     case 'sqrtscale'
%                         obj(o).F = obj(o).F ./ sqrt(obj(o).f_c) ./ sqrt(siz(1)*siz(2))
%                     case 'none'
%                     otherwise
%                         error('OrientationSpaceFilter:setupFilterNormEnergy', ...
%                             'Invalid normEnergy property');
%                 end
%             end
        end
    end
    
end

function radialKernel = radialWaveletKernel(f_c, coords)
    f_c = shiftdim(f_c(:),-2) / 0.25;
    radialKernel = pyramid.nopi.raisedCosLogFreqFilterBandpass(bsxfun(@rdivide,coords.f,f_c));
end