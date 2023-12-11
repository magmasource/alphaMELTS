classdef MELTSevent < event.EventData
    
    properties
        funcName
    end
    
    methods
        function data = MELTSevent(func)
            data.funcName = func;
        end
    end
    
end
