classdef Measure < handle
    properties (SetAccess = protected)
        mu = [];    % weights of the measure  
        coord = []; % d-dimensional coordinates of the support points
    end % properties
    
    properties (Dependent = true)
        n;  % number of nodes
        s;  % size of mu
        r;  % maximal radius of the nodes
    end % properties
    
    methods
      % Constructor
      function h = Measure(mu,coord)
        h.mu = mu;
        h.coord = coord;
        if sum(size(coord(:,1)) ~= size(mu))
           error('mu and xi must have same size')
        end
      end % function
      
      % Set functions
      function set.mu(h,mu)
        assert(sum(mu(:) < 0) == 0, 'mu must be nonnegative')
        assert(sum(mu(:))-1 < 1e-6, 'mu must have sum 1')
        h.mu = mu;
      end % function

      % Get functions
      function n = get.n(h)
        n = numel(h.mu);
      end % function
      
      function s = get.s(h)
        s = size(h.mu);
      end % function
      
      function r = get.r(h)
          inner_r = 0;
          for cdim=1:size(h.coord,2)
            inner_r = inner_r+ h.coord(:,cdim).^2;
          end
          r=sqrt(max(inner_r));
      end % function
    end % methods
end % classdef
