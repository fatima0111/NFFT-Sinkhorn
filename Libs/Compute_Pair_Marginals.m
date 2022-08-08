
function [marg2D] = Compute_Pair_Marginals(mus, phis, eta, output, direct_sum)
  % Input:
  % mus             Instances of class Measure
  % phis            Initializers of iteration
  % eta             Entropy regularization parameter
  % output          Particles position at final time output = sigma(input)
  % direct_sum      if false then use fast-sum otherwise use direct method
 
  for k=1:size(phis)
    if sum(size(mus{1}.coord(:,1)) ~= size(mus{1}.mu)) 
            error('mu and xi must have same size')
    end
  end
  
 % Number of marginals
  K= size(phis,2);

  % Parameters for fastsum
  d = size(mus{1}.coord,2);  % space dimension
  kernel = 'gaussian';  % kernel  K(x) = EXP(-x^2/c^2) 
  c = sqrt(eta);        % kernel parameter
  p = 3;                % degree of smoothness of regularization
  flags = 0;            % flags (could be EXACT_NEARFIELD or NEARFIELD_BOXES)
  n = 156;              % expansion degree
  if eta < 5e-3 %Truncate                    
      eps_I = 3/400; % inner boundary = 0 because kernel is continuous around 0
      eps_B = 0;
      flags = EXACT_NEARFIELD + NEARFIELD_BOXES;

  elseif eta >= 5e-3 && eta <= 0.05 % No truncation No regularization
       eps_I = 0;     % inner boundary = 0 because kernel is continuous around 0
       eps_B =0;      % outer boundary

  else % regularization at the boundary  
     eps_I = 0; %0.03;    % inner boundary = 0 because kernel is continuous around 0
     eps_B = 3/32;   % outer boundary
  end
  m = p;              % cut-off parameter for internal NFFT
  nn_oversampled=2*n; % oversampling factor for internal NFFT
  
  tic;
  % Rescale nodes because fastsum requires them in ball of radius
  % 1/4-eps_B/2
  inner_r =0;
  for cdim=1:size(output,2)
       inner_r = inner_r+ output(:,cdim).^2;
  end
  r_max = sqrt(max(inner_r));
  for k=1:K
    if mus{k}.r>r_max
        r_max = mus{k}.r;
    end
  end
  scale = 2^(floor(log2((1/2-eps_B)/2*r_max)));
  
  % Initialize fastsum instances for the different transforms
  kers = cell(1, K);
  inv_kers = cell(1, K);
  for k =1:K
      kers{k} = fastsum(d,kernel,c*scale,flags,n,p,eps_I,eps_B,nn_oversampled,m);
      inv_kers{k} = fastsum(d,kernel,c*scale,flags,n,p,eps_I,eps_B,nn_oversampled,m);
    
      if k~= K-1
          if k==K
           
            kers{k}.x = output.*scale;
            inv_kers{k}.y = output.*scale;
          else
            
            kers{k}.x = mus{mod(k+1, K)}.coord.*scale;
            inv_kers{k}.y = mus{mod(k+1, K)}.coord.*scale;
          end
      else
      
         kers{k}.x = mus{K}.coord.*scale; 
         inv_kers{k}.y = mus{K}.coord.*scale;
      end
      
      kers{k}.y = mus{k}.coord.*scale;
      inv_kers{k}.x = mus{k}.coord.*scale;
  end

% Initialize Betas' and gammas' array
 betas = cell(1, K-1);
 gammas = cell(1,K-1);
 % Compute initial betak for all k=K,K-1,...,2
 for n_1 = 1:mus{1}.n
      for k=K:-1:2
          if k == K
              cbeta = exp(-sum((output(n_1,:)-mus{K}.coord).^2,2)/eta);%
              %exp(-pdist2(mus{K}.coord,output).^2/eta);
          else
              phik = phis{k+1};
               kers{k}.alpha = phik(:).*betas{k}(:,n_1);
               if direct_sum ==false
                   fastsum_trafo(kers{k}); 
                  
               else
                    fastsum_trafo_direct(kers{k});
                   
               end
    
            cbeta = real(kers{k}.f);
           cbeta(cbeta<0) = 0;
          end
         
          betas{k-1}(:,n_1) = cbeta;  
      end 
 end


  % Compute first marginal divided pointwise with dual vector u1
  for n_1 = 1:mus{1}.n
    var = phis{2}(:).*betas{1}(:,n_1);
    Marg1 = var'*exp(-sum((mus{2}.coord-mus{1}.coord(n_1,:)).^2,2)/eta);
    phis{1}(n_1) = mus{1}.mu(n_1)./Marg1;
  end
    
  for k=2:K
      %Compute gamma_k for k=2,3, ..., K
      if k==2
            cgamma = exp(-pdist2(mus{2}.coord,mus{1}.coord).^2/eta);
      else
         gammak = zeros(mus{k}.n,mus{1}.n); 
         for i = 1:mus{1}.n
           inv_kers{k-1}.alpha = phis{k-1}(:).*gammas{k-2}(i,:)';
           if direct_sum ==false
               fastsum_trafo(inv_kers{k-1}); 
           else
               fastsum_trafo_direct(inv_kers{k-1}); 
           end
           gammak(:,i)= real(inv_kers{k-1}.f);
           gammak(gammak<0) = 0;
         end
         cgamma = gammak;
      end
      gammas{k-1} = cgamma';
      Margk = sum(betas{k-1}.*(phis{1}.*gammas{k-1})', 2);
      
      % Update the dual vector 
      phis{k}= mus{k}.mu(:)./Margk;
  end
  marg2D = cell(1,K-1);
  for k=2:K
    marg2D{k-1} = (phis{1}.*gammas{k-1})'.*(phis{k}.*betas{k-1});
 end
  
time = toc; 
fprintf('Computation time pair-marginal plan   %g\n', time);

end % function