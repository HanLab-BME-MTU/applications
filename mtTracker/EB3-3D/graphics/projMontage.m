function three=projMontage(XY,XZ,YZ)
  stripeSize=4;
  threeTop = [XY, zeros(size(XY,1), stripeSize,3), zeros(size(XY,1), size(YZ,2),3)];
  threeBottom = [XZ, 0*ones(size(XZ,1),+stripeSize,3),YZ];
  three = [threeTop; ones(stripeSize, size(XY,2)+size(YZ,2)+stripeSize,3); threeBottom];
