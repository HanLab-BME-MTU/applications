function newSurf = createAndSetSurfaces(this,surfaces,normals,timeIndices,name,color)

nSurf = numel(surfaces);

newSurf = this.mImarisApplication.GetFactory.CreateSurfaces;

newSurf.SetName(name)
newSurf.SetColorRGBA(this.mapRgbaVectorToScalar(color));

if isempty(normals) && isfield(surfaces,'normals')
    normals = arrayfun(@(x)(x.normals),surfaces,'Unif',false);
else
    error('You must input surface normals!')
end

if nSurf == 1
    newSurf.AddSurface(surfaces.vertices,surfaces.faces-1,normals,timeIndices)
else
    nVertPer = arrayfun(@(x)(size(x.vertices,1)),surfaces);
    nTriPer = arrayfun(@(x)(size(x.faces,1)),surfaces);
    allVert = vertcat(surfaces(:).vertices);
    allTri = vertcat(surfaces(:).faces)-1;
    allNorm = vertcat(normals{:});
    
    newSurf.AddSurfacesList(allVert,nVertPer,allTri,nTriPer,allNorm,timeIndices)
end

this.mImarisApplication.GetSurpassScene.AddChild(newSurf,-1);