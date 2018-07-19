function pdValues = LG_rag_getPdValues(pdHandles)
%LG_rag_getPdValues reads the PD-values from their handles

pdValues = get(pdHandles,'Value');
pdValues = reshape([pdValues{:}],size(pdHandles));