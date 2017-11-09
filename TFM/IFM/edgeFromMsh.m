        function edge = edgeFromMsh(msh,numEdges)
        % The geometry and mesh, 'msh' created by 'getMeshFromMask' has four boundary edges
        % parameterized by arclength. So, first identify the four edges from msh.e.
        % We use a structure named 'edge' to store all these mesh structure info and
        % displacement. It has the following fields:
        % 'vertEI' : Index into 'msh.e(3,:)' whose arclength is 0. It identifies the
        %            vertex of the edge.
        % 'I'      : Index into 'msh.p' to get the real coordinates of boundary points.
        % 'bndEI'  : The index of boundary elements that belong to one edge.
        % 'endI    : The index of the boundary element that is at the end of one
        %            edge.
        % 'arcLen' : The arclength of each edge.
        % 'bndP'   : The coordinates of boundary points on each edge.
        % 'bndS'   : The arclength parameters of the boundary points on each edge.
        % 'ppX'    : The pp-form of the spline interpolation of X-coordinate of 'bndP'.
        % 'ppY'    : The pp-form of the spline interpolation of Y-coordinate of 'bndP'.
        % 'dispV'  : The coordinates of the base and the point end of the displacement
        %            vectors on each edge in the formate [x0 y0 x1 y1].
        % 'U1'     : The first and 
        % 'U2'     : the second components of the displacement vectors on each edge.
        % 'UC1'    : For debugging. Has the same structure as 'edgeU1'.
        % 'UC2'    : For debugging. Has the same structure as 'edgeU2'.
        % 'ppU1'   : The pp-form of the spline interpolation of 'U1'.
        % 'ppU2'   : The pp-form of the spline interpolation of 'U2'.
        edge = struct('bndEI',{},'bndS',{},'endI',{},'arcLen',{},'I',{},'bndP',{},'ppX',{},'ppY',{});
        for k = 1:numEdges
         %Extract the arclenth parameters.
         edge(k).bndEI = find(msh.e(5,:)==k);
         edge(k).bndS  = msh.e(3,edge(k).bndEI);

         [edge(k).bndS,sortI] = sort(edge(k).bndS);

         edge(k).endI   = edge(k).bndEI(sortI(end));
         edge(k).bndS   = [edge(k).bndS msh.e(4,edge(k).endI)];
         edge(k).arcLen = edge(k).bndS(end);

         %Get the index of the boundary points.
         edge(k).I = msh.e(1,edge(k).bndEI);
         edge(k).I = [edge(k).I(sortI) msh.e(2,edge(k).endI)];

         edge(k).bndP  = [msh.p(1,edge(k).I); msh.p(2,edge(k).I)];
         bndSKnt = augknt(edge(k).bndS,2);
         edge(k).ppX = spapi(bndSKnt,edge(k).bndS, edge(k).bndP(1,:));
         edge(k).ppY = spapi(bndSKnt,edge(k).bndS, edge(k).bndP(2,:));
         %edge(k).ppX = spline(edge(k).bndS, edge(k).bndP(1,:));
         %edge(k).ppY = spline(edge(k).bndS, edge(k).bndP(2,:));
        end
