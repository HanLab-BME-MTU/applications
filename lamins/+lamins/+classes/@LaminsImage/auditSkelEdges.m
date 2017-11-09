function [ A ] = auditSkelEdges( skel, I )
%auditSkelEdges Check to see if skeleton edges are reasonable

endpts = bwmorph(skel,'endpoints');
branchpts = bwmorph(skel,'branchpoints');
branchpt.idx = find(branchpts);
endpt.idx = find(endpts);

vertices = branchpts | endpts;
[edges_cc, vertices_cc, pairs] = bwtrace(skel,vertices);

edges_rp = regionprops(edges_cc,I,'MaxIntensity','MeanIntensity','MinIntensity');

A.rp = edges_rp;
A.cc = edges_cc;
A.width = ([A.rp.MaxIntensity]-[A.rp.MinIntensity])./[A.rp.MeanIntensity];

thresh = thresholdRosin(A.width);
filter = (A.width < thresh | [A.rp.MeanIntensity] > 0.5) & [A.rp.MinIntensity] > 0.2;
A.fcc = filtercc(A.cc,filter);
A.frp = regionprops(A.fcc,I,'MaxIntensity','MeanIntensity','MinIntensity','Area');

vert_neigh = imdilate(vertices,strel('square',3)) & skel;

bw = labelmatrix(A.fcc) > 0 | vert_neigh;
bw = bwmorph(bw,'spur',Inf);
A.bw = bw;
%A.overlay = showSkelOnIntensity(I,bw);

end

