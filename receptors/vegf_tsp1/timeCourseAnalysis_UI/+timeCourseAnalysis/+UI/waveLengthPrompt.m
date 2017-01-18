function [ emissionWL, emissionStr, canceled ] = waveLengthPrompt( )
%waveLengthPrompt prompt for wavelength selection
%
% OUTPUT
% emissionwL: Number indicating the wavelength selected or empty if canceled

% Tae Kim, originally part of UITimeCourseAnalysis
% Adapted into a function by Mark Kittisopikul, Nov 2015

    canceled = false;
    emissionWL = -1;
    emissionStr = [];
    
    emissionStr = inputdlg('Enter emission wavelength:','Enter emission wavelength:');
    if(isempty(emissionStr))
        canceled = true;
    else
        emissionWL  = str2double(emissionStr);
    end
    
    return;

    stringListWL = { ...
        '525nm : Alexa 488' ,  525
        '530nm : GFP'       ,  530 
        '590nm : Rhod Red X',  590 
        '668nm : Alexa 640' ,  668 
        '669nm : Atto 547N' ,  559 
        'Brightfield'       ,   []
        };
    userChoiceWL = listdlg('PromptString','Select wavelength:', 'SelectionMode','single', 'ListString', stringListWL(:,1));
    
    if(isempty(userChoiceWL))
        % Selection canceled
%         emissionWL = [];
        canceled = true;
    else
        % Selection made;
        emissionWL = stringListWL{userChoiceWL,2};
        if(nargout > 1)
            emissionStr = stringListWL{userChoiceWL,1};
        end
    end
%     if userChoiceWL == 1
%         emissionWL = 525;
%     end
%     if userChoiceWL == 2
%         emissionWL = 530;
%     end
%     if userChoiceWL == 3
%         emissionWL = 590;
%     end
%     if userChoiceWL == 4
%         emissionWL = 668;
%     end
%     if userChoiceWL == 5
%         emissionWL = 669;
%     end
%     if userChoiceWL == 6
%         emissionWL = [];
%         param.imageType_ = '';
%     end
    
end

