classdef PriorityQueue < handle

    properties (SetAccess = private)               
       numElements; 
       priorityList;
       valueList;
       flagMaxPriorityQueue;
    end

    methods (Access = public)       

        function obj = PriorityQueue( flagMaxPriorityQueue )

            if ~exist( 'flagMaxPriorityQueue', 'var' )
                flagMaxPriorityQueue = true;
            else
                if ~(isscalar(flagMaxPriorityQueue) && islogical(flagMaxPriorityQueue))
                    error( 'ERROR: invalid flagMaxPriorityQueue argument' );
                end
            end

            obj.flagMaxPriorityQueue = flagMaxPriorityQueue;
            obj.numElements = 0;
            obj.priorityList = {};
            obj.valueList = {};

        end

        function insert(obj, priority, value)     

            % increase the size of the array if full
            if obj.numElements > 0 && obj.numElements + 1 > numel( obj.priorityList )                                

                % double the size of the array and copy stuff
                obj.priorityList = cat(1, obj.priorityList, cell(obj.numElements, 1));
                obj.valueList = cat(1, obj.valueList, cell(obj.numElements, 1));

            end

            obj.numElements = obj.numElements + 1;

            obj.priorityList{ obj.numElements } = priority;
            obj.valueList{ obj.numElements } = value;

            obj.swim(obj.numElements);

        end

        function [priority, value] = pop( obj )

            if obj.isEmpty()
                error( 'called pop() on an empty priority queue' );
            end          

            priority = obj.priorityList{1};
            value = obj.valueList{1};

            obj.exch(1, obj.numElements);            
            obj.numElements = obj.numElements - 1;            
            obj.sink(1);

            obj.priorityList{ obj.numElements + 1 } = [];
            obj.valueList{ obj.numElements + 1 } = [];

            % halve the size of the arrays if they get one-quarter full
            if obj.numElements > 0 && obj.numElements == floor( numel( obj.priorityList ) / 4 )                

                obj.priorityList( 2 * obj.numElements + 1 : end ) = [];
                obj.valueList( 2 * obj.numElements + 1 : end ) = [];

            end

        end

        function [flagEmpty] = isEmpty( obj )        

            flagEmpty = (obj.numElements == 0);

        end

        function [qSize] = size( obj )

            qSize = obj.numElements;

        end

        function [priority, value] = peek( obj )

            if obj.isEmpty()
                error( 'requested max() of an empty priority queue' );
            end          

            priority = obj.priorityList{1};
            value = obj.valueList{1};

        end

    end    

    methods (Access = private)

        function swim(obj, elPos)

            while elPos > 1 && obj.compare(floor(elPos / 2), elPos)

                obj.exch(floor(elPos / 2), elPos);
                elPos = floor(elPos / 2);

            end

        end

        function sink(obj, elPos)

            while 2 * elPos <= obj.numElements

                j = 2 * elPos;

                if j < obj.numElements && obj.compare(j, j+1)
                    j = j + 1;
                end

                if ~obj.compare(elPos, j)
                    break;
                end

                obj.exch(elPos, j);
                elPos = j;

            end

        end

        function [blnCmpResult] = compare(obj, e1, e2)

            if obj.flagMaxPriorityQueue
                blnCmpResult = (obj.priorityList{e1} < obj.priorityList{e2});
            else
                blnCmpResult = (obj.priorityList{e1} > obj.priorityList{e2});
            end            

        end

        function exch(obj, e1, e2 )

            temp = obj.priorityList{e1};
            obj.priorityList{e1} = obj.priorityList{e2};
            obj.priorityList{e2} = temp;            

            temp = obj.valueList{e1};
            obj.valueList{e1} = obj.valueList{e2};
            obj.valueList{e2} = temp;            

        end

    end

end %classdef
