function T = direct_align(src, dst, levels, iterations, mv_type, dbg_flag, init_guess, mask)

%% ================================================================================
% direct_align 
% -------------------------------------------------
%
% Inputs:
%  ~~~~~~
% src
% dst
% levels
% iterations
% mv_type
% dbg_flag
% init_guess
% mask
%
% Outputs:
% ~~~~~~~~
% T
%
%% ================================================================================

if(~exist('init_guess','var') || isempty(init_guess))
    init_guess = eye(3);
end
T = init_guess;

if(~exist('mask','var'))
    mask = ones(size(src,1), size(src,2));
end

src = double(src);
dst = double(dst);

srcPyramid = createPyramid(src, levels);
dstPyramid = createPyramid(dst, levels);
maskPyramid = createPyramid(mask, levels);

S = eye(3); S([1,5]) = 0.5;

% Chenge the scale of the initial transform to the highest level
for l=2:levels
    T = S*T*inv(S);
end

for l=levels:-1:1
    
    T = SingleLevelFindMotion(srcPyramid{l}, dstPyramid{l}, maskPyramid{l}, T, iterations, mv_type,dbg_flag);
    if(any(isnan(T(:))) )
        return;
    end
    if(l > 1)
        T = inv(S)*T*S;      
    end
end

end % of function: direct_align

%%
function  T = SingleLevelFindMotion(src, dst, mask, T, iterations, mv_type,dbg_flag)

    center = size(mask)/2;
    changeCoordinatesTransform = eye(3);
    changeCoordinatesTransform(1:2,3) = center(1:2);
    ccti = inv(changeCoordinatesTransform);

    nEq = numberOfEquation(mv_type);

    for i=1:iterations
        inv_T = inv(T);
        if(any(isnan(inv_T(:))) || any(isinf(inv_T(:))))
            %make sure T is a valid transform.exit if not
            return;
        end
        tempMask = wrap_affine(mask, inv_T);

        DX = filter2([-0.5,0,0.5],src);
        DY = filter2([-0.5,0,0.5]',src);
        MD = wrap_affine(dst, inv_T) - src;

        [mat, vect] = FindMotion(DX, DY, MD,tempMask, mv_type, ccti);

        R = mat\vect;
        T0 = eye(3);
        if(strcmp(mv_type,'trans'))
            T0(1,3) = R(1);
            T0(2,3) = R(2);
        elseif(strcmp(mv_type,'affine'))
            T0(1,1) = R(1);
            T0(1,2) = R(2);
            T0(1,3) = R(3);
            T0(2,1) = R(4);
            T0(2,2) = R(5);
            T0(2,3) = R(6);
        end
        T0 = changeCoordinatesTransform*T0*inv(changeCoordinatesTransform);
        T = T*T0;

        if(any(isnan(T(:))))
            %make sure T is a valid transform.exit if not
            return;
        end

        if dbg_flag
            imshow(MD/255+0.5);
            pause(0.1);
        end

        disp = TransformMaxDisparity(T0, size(src));
        if(disp < 0.001)
            break;
        end

    end

end % of function: SingleLevelFindMotion

%%
function disparity = TransformMaxDisparity(T, vid_size)
s = vid_size;
p = ones(3,4);
i=0:3;
p(1,mod(i,2)==1) = s(1);
p(2,mod( floor(i/2),2)==1) = s(2);

q = T*p;

disparity = max( sum((p-q).^2));

end % of function: TransformMaxDisparity

%%
function [mat,vect] = FindMotion(DX, DY, MD, mask, mv_type, ccti)

mask([1,2,end-1,end],:)=0;
mask(:,[1,2,end-1,end])=0;

[X,Y] = meshgrid(1:size(DX,2), 1:size(DX,1));
X = X(:) + ccti(1,3);
Y = Y(:) + ccti(2,3);

mask = mask(:);
DX = DX(:).*mask;
DY = DY(:).*mask;
MD = MD(:).*mask;

if(strcmp(mv_type,'trans'))
    A = [DX,DY];
    nit = -MD;
elseif(strcmp(mv_type,'affine'))
    A = [X.*DX, Y.*DX, DX, X.*DY, Y.*DY, DY];
    nit = -MD + X.*DX + Y.*DY;
end

mat = A'*A;
vect = A'*nit;

end % of function: FindMotion

%%
function nEq = numberOfEquation(mv_type)

if(strcmp(mv_type,'trans')),nEq=2;end
if(strcmp(mv_type,'affine')),nEq=6;end

end % of function: numberOfEquation

%%
function pyr = createPyramid(im, levels)

pyr = cell(levels,1);
pyr{1} = im;

for i=2:levels
    pyr{i} = reduce(pyr{i-1});
end

end % of function: createPyramid

%%
function  out=reduce(image)
    %--------------------------------------------------------------------
    %   Reduce the image size to half by blurring  & sub-sampling
    %
    %   image    : the input image
    %
    %  function out = Reduce(image)
    %--------------------------------------------------------------------

    % the blurring gaussian kernel
    % [c b a b c]
    % b = 1/4
    % a = 0.4
    % c = b - a/2
    kernel = [0.05 0.25 0.4 0.25 0.05];
    kernel = kernel'*kernel;

    % filter
    B = filter2 (kernel, image, 'valid');
    % sub-sampling the image
    out = B(1:2:end,1:2:end);

end % of function: Reduce
