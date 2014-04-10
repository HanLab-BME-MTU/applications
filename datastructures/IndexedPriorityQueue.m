classdef IndexedPriorityQueue < handle
    
   properties (SetAccess = private) 
        numElements;
        priorityList;
        indexList;    
        indexPosList;
        flagMaxPriorityQueue;
   end
   
   methods (Access = public)
       
       function obj = IndexedPriorityQueue( maxIndex, flagMaxPriorityQueue )
           
            if ~exist( 'flagMaxPriorityQueue', 'var' )
                flagMaxPriorityQueue = true;
            else
                if ~(isscalar(flagMaxPriorityQueue) && islogical(flagMaxPriorityQueue))
                    error( 'ERROR: invalid flagMaxPriorityQueue argument' );
                end
            end
           
            obj.flagMaxPriorityQueue = flagMaxPriorityQueue;
            obj.numElements = 0;
            
            obj.priorityList = [];
            obj.indexList = [];            
            obj.indexPosList = zeros(1, maxIndex);
           
       end
       
        function insert(obj, index, priority)     

            % make sure if index is not present already
            if obj.contains(index)
                error( 'ERROR: index is already present in the priority queue' );
            end
            
            % increase the size of the array if full
            if obj.numElements > 0 && obj.numElements + 1 > numel( obj.priorityList )                                

                % double the size of the array and copy stuff
                obj.priorityList = cat(1, obj.priorityList, zeros(obj.numElements, 1));
                obj.indexList = cat(1, obj.indexList, zeros(obj.numElements, 1));

            end
            
            obj.numElements = obj.numElements + 1;

            obj.priorityList( obj.numElements ) = priority;
            obj.indexList( obj.numElements ) = index;   
            obj.indexPosList( index ) = obj.numElements;
            
            obj.swim(obj.numElements);

        end

        function [index, priority] = pop( obj )

            if obj.isEmpty()
                error( 'called pop() on an empty priority queue' );
            end          

            priority = obj.priorityList(1);
            index = obj.indexList(1);

            obj.exch(1, obj.numElements);            
            obj.numElements = obj.numElements - 1;            
            obj.sink(1);

            obj.priorityList( obj.numElements + 1 ) = [];
            obj.indexList( obj.numElements + 1 ) = [];
            obj.indexPosList( index ) = 0;

            % halve the size of the arrays if they get one-quarter full
            if obj.numElements > 0 && obj.numElements == floor( numel( obj.priorityList ) / 4 )                

                obj.priorityList( 2 * obj.numElements + 1 : end ) = [];
                obj.indexList( 2 * obj.numElements + 1 : end ) = [];

            end

        end

        function [flagEmpty] = isEmpty( obj )        

            flagEmpty = (obj.numElements == 0);

        end

        function [qSize] = size( obj )

            qSize = obj.numElements;

        end

        function [index, priority] = peek( obj )

            if obj.isEmpty()
                error( 'requested max() of an empty priority queue' );
            end          

            priority = obj.priorityList(1);
            index = obj.indexList(1);

        end
        
        function [flagPresent] = contains( obj, index )
            
            flagPresent = obj.indexPosList( index ) > 0;
            
        end
        
        function update(obj, index, newPriority)
            
            if ~contains( obj, index )
               error( 'ERROR: attempted to update priority of an index which is not present in the priority queue' );
            end
            
            indexPos = obj.indexPosList( index );
            curPriority = obj.priorityList( indexPos );
            
            if obj.comparePriorities(newPriority, curPriority)

               obj.priorityList( indexPos ) = newPriority;
               obj.sink( indexPos );
               
            else

               obj.priorityList( indexPos ) = newPriority;
               obj.swim( indexPos );
                
            end
            
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
        
        function [blnCmpResult] = comparePriorities(obj, p1, p2)

            if obj.flagMaxPriorityQueue
                blnCmpResult = (p1 < p2);
            else
                blnCmpResult = (p1 > p2);
            end            
            
            
        end
        
        function [blnCmpResult] = compare(obj, e1, e2)

            blnCmpResult = obj.comparePriorities(obj.priorityList(e1), obj.priorityList(e2)); 

        end

        function exch(obj, e1, e2 )

            e1_index = obj.indexList(e1);
            e1_priority = obj.priorityList(e1);

            e2_index = obj.indexList(e2);
            e2_priority = obj.priorityList(e2);
            
            obj.priorityList(e1) = e2_priority;
            obj.priorityList(e2) = e1_priority;            

            obj.indexList(e1) = e2_index;
            obj.indexList(e2) = e1_index;            

            obj.indexPosList( e1_index ) = e2;
            obj.indexPosList( e2_index ) = e1;
            
        end           
           
   end
    
end