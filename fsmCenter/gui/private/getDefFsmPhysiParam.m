function physiParam = getDefFsmPhysiParam

physiParam.NA            = 1.4;
physiParam.waveLen       = 600; %Unit, nm.
physiParam.bitDepth      = 14;
physiParam.pixelSize     = 67; %Unit, nm.
physiParam.frameInterval = 10; %Unit, sec.
physiParam.psfSigma      = round(0.21*physiParam.waveLen/physiParam.NA ...
                           /physiParam.pixelSize*1000)/1000;

