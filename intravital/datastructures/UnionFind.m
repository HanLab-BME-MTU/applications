classdef UnionFind < handle
% The Weighted Quick Union implementation of the union-find data structure with path compression 

    methods (Access = public)
    
        % create and intialize a union-find data structure
        function this = UnionFind( numObjects )
            
           this.objVec = 1:numObjects; 
           this.subtreeSizeVec = ones(1, numObjects);
           
        end
        
		% finds and returns the label of the connnected component of <obid>
        function compid = find(this, obid)            
            
            compid = obid;
            while this.objVec( compid ) ~= compid                
                compid = this.objVec( compid );
            end
        
            % path compression -- flatten the tree
            curObid = obid;
            while this.objVec( curObid ) ~= curObid
                pid = this.objVec( curObid );                
                this.objVec( curObid ) = compid;
                curObid = pid;
            end
            
        end
        
		% joins the connnected component of <obid1> and <obid2> 
        function union(this, obid1, obid2)
                        
            root1 = this.find( obid1 );
            root2 = this.find( obid2 );
            
            if this.subtreeSizeVec(root1) > this.subtreeSizeVec(root2)               
                this.objVec( root2 ) = root1;
                this.subtreeSizeVec(root1) = this.subtreeSizeVec(root1) + this.subtreeSizeVec(root2);
            else                
                this.objVec( root1 ) = root2;
                this.subtreeSizeVec(root2) = this.subtreeSizeVec(root2) + this.subtreeSizeVec(root1);
            end
            
        end
        
		% checks if <obid1> and <obid2> reside in the same connected component
        function [blnAreConnected] = connected(this, obid1, obid2)
            
            blnAreConnected = this.find( obid1 ) == this.find( obid2 );
            
        end
        
    end
    
    properties (SetAccess = private)        
        objVec;
        subtreeSizeVec;
    end
    
end % classdef