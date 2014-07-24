function plotAllPropertiesRelToActivityOnset(sptPropInWindow)

% numType = length(sptPropInWindow);

for iType = 1 : 1
    
   sptPropInWindowType = sptPropInWindow(iType);
   
   globalField = fieldnames(sptPropInWindowType);
   numGlobalField = length(globalField);
   
   for iGlobalField = 1 : numGlobalField
       
       propertiesGlobalField = sptPropInWindowType.(globalField{iGlobalField});
       
       localField = fieldnames(propertiesGlobalField);
       numLocalField = length(localField);
       
       for iLocalField = 1 : numLocalField
           
           propertyLocalField = propertiesGlobalField.(localField{iLocalField});
           
           figName = ['Prot. Type ' num2str(iType) ' ' globalField{iGlobalField} ' ' localField{iLocalField}];
           plotSptRelToActivityOnsetAdaptiveWindows(propertyLocalField,[],[],figName)
           
       end
       
   end
    
end