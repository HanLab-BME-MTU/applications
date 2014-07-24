function nbList=constructNeighborList(grid,varargin)
%CONSTRUCTNEIGHBORLIST returns neighboring bins of a regular lattice
%
%   required input arguments:
%       grid -> grid dimensions
%   
%   optional input arguments:
%       mode -> mode can be 'periodic','truncated', default 'periodic'
%
%   output:
%       nbList -> list of neighbors for each bin in grid
%
%   US, 2012/11/20
%

ip=inputParser;

ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('grid',@isnumeric);
ip.addOptional('mode', 'periodic', @(x) any(strcmpi(x, {'periodic', 'truncated'})));

ip.parse(grid,varargin{:});

grid=ip.Results.grid;
mode=ip.Results.mode;

nbList=cell(grid);

for n=1:prod(grid)
    idc=1;
    clique=NaN(8,2);

    [idx,idy]=ind2sub(grid,n);
    % row above current bin: idx -> idx-1
    xx=getNeighbor(idx,-1,grid(1),'mode',mode);
    for k=-1:1
        yy=getNeighbor(idy,k,grid(2),'mode',mode);
        clique(idc,:)=[xx,yy];
        idc=idc+1;
    end
    
    % row below current bin: idx -> idx+1
    xx=getNeighbor(idx,1,grid(1),'mode',mode);
    for k=-1:1
        yy=getNeighbor(idy,k,grid(2),'mode',mode);
        clique(idc,:)=[xx,yy];
        idc=idc+1;
    end
    % bins left and right of current bin
    xx=idx;
    yy=getNeighbor(idy,-1,grid(2),'mode',mode);
    clique(idc,:)=[xx,yy];
    idc=idc+1;
    yy=getNeighbor(idy,1,grid(2),'mode',mode);
    clique(idc,:)=[xx,yy];
    % filter invalid neighbors
    idc=true(8,1);
    for k=1:8
        idc(k)=~any(clique(k,:) == -1);
    end
    % store in nbList
    nbList{idx,idy}=clique(idc,:);   
end


end