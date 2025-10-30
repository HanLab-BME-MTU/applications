function mustBeField(a)
    if ~isfield(a,'pos') || ~isfield(a,'vec')
        eidType = 'mustBeField:notField';
        msgType = 'Input field is incompatible, please input a structure which contains the pos and vec fields, and only those fields.';
        error(eidType,msgType)
    end
end