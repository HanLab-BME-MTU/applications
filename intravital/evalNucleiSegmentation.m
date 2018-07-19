function [ cmdOutput, imLabelCellSeg, imCellSeedPoints ] = evalNucleiSegmentation(segFunc, im)

    [cmdOutput, imLabelCellSeg, imCellSeedPoints] = evalc( 'segFunc(im)' );

end