function tracksNA = getTracksNAFromMD(MD)

iAdhProc = MD.getProcessIndex('AdhesionAnalysisProcess');
adhAnalProc = MD.getProcess(iAdhProc);
p = adhAnalProc.funParams_;
tracksNA=adhAnalProc.loadChannelOutput(p.ChannelIndex,'output','tracksNA');

