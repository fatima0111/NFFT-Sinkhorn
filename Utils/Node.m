classdef Node < handle
    properties (SetAccess = protected)
        key = -1;                  % Node Id
        marginal = Measure.empty;  % Node as discrete probability measure
    end % properties
    
    properties
        parent = Node.empty;    % Parent of the node  
        children = [];          % Collection of Children 
        is_leaf = false;        % Boolean to specify wether the node is a leaf 
        barycenter_weights = 1; % weights of the Node in the generalized barycenter problem
    end % properties

    
    methods
      % Constructor
      function N = Node(key, mu, coord)
        N.key = key;
        N.marginal = Measure(mu, coord);
        if key <1 
           error('the key must be larger than 1')
        end
      end % end function
      
      % Set functions
      function set.key(N,key)
        N.key = key;
      end % end function
      
      function set.children(N,children)
        N.children = children;
      end % end function
      
      function set.is_leaf(N,Is_leaf)
        N.is_leaf = Is_leaf;
      end % end function
   
      function set.barycenter_weights(N,b_weights)
        N.barycenter_weights = b_weights;
      end % end function
      
      % Get function
      function parent= get.parent(N)
          parent = N.parent;
      end

      function children= get.children(N)
          children = N.children;
      end % end function

    end % methods
end % classdef
