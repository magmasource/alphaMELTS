classdef MELTSmap < containers.Map
    %MELTSmap Subclass of containers.Map that takes string arrays as KEYS and allows 2-D indexing.
    %   Detailed explanation goes here    
    methods
        function obj = MELTSmap(varargin)
            %MELTSmap Construct an instance of this class with 'KeyType' = string array, and 'ValueType' = 'any'
            %   
            %   mapObj = MELTSmap constructs an empty MELTSmap container mapObj
            %
            %   mapObj = MELTSmap(keySet,valueSet) constructs a Map that contains one or more values and a unique key for each value.

            obj@containers.Map();
            if nargin
                obj.subsasgn(struct('type', '()', 'subs', varargin{1}), varargin{2});
            end
            
        end
        
        function obj = subsasgn(obj, keys, values)
            %subsasgn Subscripted assignment.
            %   Detailed explanation goes here
            
            keys = keys.subs;
            if iscell(keys) && isstring(keys{:})
               keys = keys{:};
            end
            if isstring(keys) || ischar(keys)
                keys = cellstr(keys);
            end            
            
            if length(keys) > 1 && ismatrix(values)
                if numel(keys) == numel(values)
                    keys = reshape(keys, [], 1);
                    if size(keys) == size(values)
                        values = reshape(values, [], 1);
                    else
                        values = reshape(values, 1, []);
                    end
                    values = num2cell(values);
                elseif isrow(keys)
                    values = mat2cell(values, size(values, 1), ones(1, length(keys)));
                elseif iscolumn(keys)
                    values = mat2cell(values, ones(length(keys), 1), size(values, 2));
                end                
            elseif length(keys) == 1 && ~iscell(values)
                values = {values};
            end
            
            assert(length(keys) == length(values), 'The number of keys and values must be the same.')
            for i = 1:length(keys)
                subsasgn@containers.Map(obj, struct('type', '()', 'subs', keys{i}), values{i});
            end
            
        end

        function value = subsref(obj, keys)
            %subsref Subscripted reference.
            %   Detailed explanation goes here
            % Assumes all vectors the same length - should have a warning
            
            if isscalar(keys)
                if strcmp(keys.type, '()')
                    keys = keys.subs;
                    if iscell(keys) && isstring(keys{:})
                        keys = keys{:};
                    end
                    if isstring(keys) || ischar(keys)
                        keys = cellstr(keys);
                    end
                    
                    value = [];
                    for i = 1:length(keys)
                        el = subsref@containers.Map(obj, struct('type', '()', 'subs', keys{i}));
                        if isrow(keys); value = [value el];
                        else; value = [value; el]; end
                    end
                else
                    value = subsref@containers.Map(obj, keys);
                    if iscellstr(value); value = string(value); end
                end
            else
                value = subsref@containers.Map(obj, keys(1));
                value = subsref(value, keys(2));
                if ischar(value); value = string(value); end
            end
        end

        
    end
end

