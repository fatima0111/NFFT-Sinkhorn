function [instance_tree, phis] =  Generate_tree(K)
    %
    % Input
    % K            The number of nodes of the tree
    keys = zeros(1,K);
    phis = cell(1, K);
    ms = cell(1,K);
    coords = cell(1, K);
    edges = cell(1, K-1);

    for k=1:K
        C = sort(unifrnd(-0.5,0.5, n_p,1));
        ms{k} = zeros(n_p,1)+ 1/n_p;
        mu = Measure(zeros(n_p,1)+ 1/n_p,C);
        keys(k) = k;
        phis{k} = ones(mu.s)/mu.n;
        coords{k} = C;
    end
    edges{1} = [1,2];
    height = floor(log2(K-1));
    c_nodes = {2};
    c_nodes_new = {};
    l=3;
    while height ~=0
       iter_l = 1;
       for n_nodes = 1:size(c_nodes,2)
           if l <=K-1
               edges{l-1} = [c_nodes{n_nodes},l];
               edges{l} = [c_nodes{n_nodes},l+1];
               c_nodes_new{iter_l} = l;
               c_nodes_new{iter_l+1} = l+1;
               l = l+2;
               iter_l = iter_l+2;
           elseif l == K 
               iter_l = iter_l +1;
               edges{l-1} = [c_nodes{n_nodes},l];
               c_nodes_new{iter_l} = l; 
               break;
           end
       end
       c_nodes = c_nodes_new;
       c_nodes_new={};
       height = height-1;
   end
   instance_tree = Tree(1,keys, edges, ms, coords);
end % end function