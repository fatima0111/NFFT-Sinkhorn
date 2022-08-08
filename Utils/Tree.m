classdef Tree < handle
    properties (SetAccess = protected)
        root  = 1;       % key of the root of the tree
        edges = [];      % set of edges of the tree
        nodes = [];      % set of nodes of the tree
        free_nodes = []; % set of node measures that are not constrained   
    end % properties
    

    methods
      % Constructor
      function t = Tree(root_key,keys, edges, mus, coords, free_contrained, b_weights)
        keys = sort(keys); 
        t.edges = edges;
        disp(t.edges);
        t.root = Node(root_key, mus{root_key}, coords{root_key});
        for e = 1:size(edges, 2) % Checks if for each p(k)<k
            if edges{e}(1)> edges{e}(2)
                error('The key of node %g is smaller than his parent %g', edges{e}(2), edges{e}(1));
            end 
        end 
        for k=keys
            nodes(k) = Node(keys(k), mus{k},coords{k});
        end
        t.nodes = nodes;
        if size(edges,2) ~= size(nodes)-1
           error('the graph is not a tree');
        end
        set_nodes_ds(t, t.root);
        iter = 0;
        iter2 = 0;
        if nargin > 5 && size(free_contrained,2)>=1

          for k=1:size(t.nodes,2)
              if ismember(k, free_contrained)
                  iter = iter +1;
                  free_nodes(iter)=t.nodes(k);
              else
                  iter2 = iter2+1;
                  t.nodes(k).barycenter_weights = b_weights(iter2);
              end
          end
          t.free_nodes= free_nodes;
        end
       
      end % end Constructor
      
      % Set recursively children and parent of node  
      function set_nodes_ds(t,node)
          indices = find(cellfun(@(edge)edge(1),t.edges)==node.key);
          if isempty(indices) ==0  % not empty
                children_keys = cellfun(@(edge)edge(2),t.edges(indices));
                t.nodes(node.key).children = t.nodes(children_keys);
                for j=1:size(children_keys,2)
                    child_key = children_keys(j);
                    t.nodes(child_key).parent = node;
                    set_nodes_ds(t, t.nodes(child_key));
                end
          else
              t.nodes(node.key).is_leaf = true;
           end
      end % end function

      % Set functions
      function set.nodes(t,nodes)
        t.nodes = nodes;
      end % end function
    end % methods
end % classdef
