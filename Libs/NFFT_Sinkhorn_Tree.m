
function [phis, ds, averaged_time] = NFFT_Sinkhorn_Tree(tree, phis, eta, direct_sum, max_it, stop_criterion)
  % Input:
  % mus             Instances of class Measure
  % phis              Initializers of iteration
  % max_it          Maximum number of iterations
  % eta             Entropy regularization parameter
  % stop_criterion  Stopping criterion for the algorithm
  
  if nargin < 6
    stop_criterion = 1e-20;
  end
  
  for k=1:size(phis)
    if sum(size(tree.nodes(k).marginal.coord(:,1)) ~= size(tree.nodes(k).marginal.mu)) 
            error('mu and xi must have same size')
    end
  end
  
  % Number of marginals
  K= size(phis,2);

  % Parameters for fastsum
  d = size(tree.nodes(1).marginal.coord,2); % space dimension
  kernel = 'gaussian';  % kernel  K(x) = EXP(-x^2/c^2) 
  c = sqrt(eta);        % kernel parameter
  p = 3;                % degree of smoothness of regularization
  flags = 0;            % flags (could be EXACT_NEARFIELD or NEARFIELD_BOXES)
  n = 156;              % expansion degree
  eps_I = 0;            % inner boundary = 0 because kernel is continuous around 0
  if eta > 0.01
    eps_B = 1/32;       % outer boundary
  else
    eps_B = 0;
  end
  m = p;              % cut-off parameter for internal NFFT
  nn_oversampled=2*n; % oversampling factor for internal NFFT
  
 
  % Rescale nodes because fastsum requires them in ball of radius
  r_max = 0; 
  for edg_ind=1:size(tree.edges)
    if (tree.nodes(tree.edges{edg_ind}(1)).marginal.r+tree.nodes(tree.edges{edg_ind}(2)).marginal.r)>r_max
        r_max =tree.nodes(tree.edges{edg_ind}(1)).marginal.r+tree.nodes(tree.edges{edg_ind}(2)).marginal.r;
    end
  end
  scale = 2^(floor(log2((1/4-eps_B)/r_max)));
     

  % Initialize fastsum instances for the different transforms
  kers = cell(1, K-1);
  inv_kers = cell(1, K-1);

  for k=1:size(tree.nodes,2)-1
      if k==1
          lambda_k = sqrt(tree.nodes(1).barycenter_weights);
      else
          lambda_k = sqrt(tree.nodes(k+1).barycenter_weights);
      end
      kers{k} = fastsum(d,kernel,(c*scale)/lambda_k,flags,n,p,eps_I,eps_B,nn_oversampled,m);
      kers{k}.x = tree.nodes(k+1).marginal.coord.*scale;
      kers{k}.y = tree.nodes(k+1).parent.marginal.coord.*scale; 
      
      inv_kers{k} = fastsum(d,kernel,(c*scale)/lambda_k,flags,n,p,eps_I,eps_B,nn_oversampled,m);
      inv_kers{k}.x = tree.nodes(k+1).parent.marginal.coord.*scale;
      inv_kers{k}.y =  tree.nodes(k+1).marginal.coord.*scale;
  end


% Initialize Betas' and gammas' array
 betas = cell(1, K-1); % the betas are vectors
 % Compute initial betak 
 for k=K:-1:2
     c_betas = ones(size(tree.nodes(k).marginal.mu));    
     for t_ind=1:size(tree.nodes(k).children,2)% If child node k is not a leaf
         t = tree.nodes(k).children(t_ind).key;
         c_betas = c_betas.*betas{t-1};
     end 
     kers{k-1}.alpha = phis{k}(:).*c_betas;
     % Matrix multiplication
     if direct_sum ==false
         fastsum_trafo(kers{k-1});
     else
         fastsum_trafo_direct(kers{k-1});
     end
     betas{k-1} = real(kers{k-1}.f);
     betas{k-1}(betas{k-1}<0)=0; %abs(betas{k-1}(betas{k-1}<0));
 end

 phis_old = phis;
 
 d = -50;
 d_old = -100;
 iteration = 0;
 %ds = zeros(max_it,1);
 norm = 1;
 % Iteratively compute maximizing sequence of dual vectors u1, u2, ..., uK
 while norm >stop_criterion  && (iteration ==0 || iteration <max_it) %abs(d-d_old)
      tic;
      d_old = d;
      phis_old = phis;
      d = 0;
      %Compute gamma_k for k=1,2, ..., K
      gammas = cell(1,K);
      Marg_root = ones(size(tree.nodes(1).marginal.mu));
      for  k=1:K
            if k==1 
                 gammas{1}=ones(size(tree.nodes(1).marginal.mu));
                 for j=1:size(tree.nodes(1).children,2)
                     t = tree.nodes(1).children(j).key;
                     Marg_root = Marg_root.*betas{t-1};
                 end
                 d=d+log(phis{1}(:)'+1e-15)*tree.nodes(1).marginal.mu(:);
                 phi1_old = phis{1};
                 if ~ismember(tree.nodes(1), tree.free_nodes)
                     phis{1} = zeros(size(tree.nodes(1).marginal.mu(:)));
                     mu_k = tree.nodes(1).marginal.mu(:);
                     phi_k  = mu_k(mu_k>0)./ Marg_root(mu_k>0);
                    phis{1}(mu_k>0) = phi_k;
                 end
                 d2 = phi1_old(:)'*Marg_root; 
            else
                k_parent = tree.nodes(k).parent;
                p_k = k_parent.key;
                c_betas = ones(size(tree.nodes(p_k).marginal.mu));
                if size(k_parent.children,2)>1
                     for j=1:size(k_parent.children,2)        
                         t = k_parent.children(j).key;
                         if t~=k
                            c_betas = c_betas.*betas{t-1};
                         end
                     end
                end
                var = gammas{p_k}.*c_betas;
                inv_kers{k-1}.alpha = phis{p_k}.*var;

                % Matrix-vector multiplication
                if direct_sum ==false
                     fastsum_trafo(inv_kers{k-1});
                else
                     fastsum_trafo_direct(inv_kers{k-1});
                end
                gammas{k} = real(inv_kers{k-1}.f);
                gammas{k}(gammas{k}<0)= 0; %abs(gammas{k}(gammas{k}<0));

                % Update dual potential phi_k
                c_betas = ones(size(tree.nodes(k).marginal.mu(:)));
                for j=1:size(tree.nodes(k).children,2)
                    t = tree.nodes(k).children(j).key;
                    c_betas = c_betas.*betas{t-1};
                end
                d = d + log(phis{k}(:)'+1e-15)*tree.nodes(k).marginal.mu(:); %Compute d before marginal potential
                if ~ismember(tree.nodes(k), tree.free_nodes)
                    %us{k} = nodes(k_ind).marginal.mu(:)./(gammas{k}.*c_betas);
                    phis{k} = zeros(size(tree.nodes(k).marginal.mu));
                    mu_k = tree.nodes(k).marginal.mu(:);
                    marg_k = (gammas{k}.*c_betas);
                    phi_k  = mu_k(mu_k>0)./marg_k(mu_k>0);
                    phis{k}(mu_k>0) = phi_k;
                end
             end  
       end 
      d = eta*(d-d2);%(u_root'*Marg_root)); 
      ds(iteration+1) = d;
      norm = 0;
      for k=1:K
          if ~ismember(tree.nodes(k), tree.free_nodes)
            norm = norm + max(abs(phis{k}-phis_old{k}));
          end
      end
      fprintf('Iteration %02d,  d = %g, norm = %g\n', iteration, d, norm);
      
      iteration = iteration+1;
      
  %Update betas 
  for k=K:-1:2
     c_betas = ones(size(tree.nodes(k).marginal.mu));    
     for t_ind=1:size(tree.nodes(k).children,2) % Go inside for-loop if k is not a leaf
         t = tree.nodes(k).children(t_ind).key;
         c_betas = c_betas.*betas{t-1};
     end 
     kers{k-1}.alpha = phis{k}(:).*c_betas;
     % Matrix multiplication
     if direct_sum ==false
         fastsum_trafo(kers{k-1});
     else
         fastsum_trafo_direct(kers{k-1});
     end
     betas{k-1} = real(kers{k-1}.f);
     betas{k-1}(betas{k-1}<0)=0; %abs(betas{k-1}(betas{k-1}<0));
 end
 
 times(iteration)=toc;
 disp(times(iteration));
 end % end while
 averaged_time = mean(times);

end % function