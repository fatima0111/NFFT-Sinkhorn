
function [phis, ds, averaged_time] = NFFT_Sinkhorn_Circle(mus, phis, eta, output, direct_sum, max_it, stop_criterion)
  % Input:
  % mus             Instances of class Measure
  % phis            Initializers of iteration
  % eta             Entropy regularization parameter
  % output          Particles position at final time output = sigma(input)
  % direct_sum      if false then use fast-sum otherwise use direct method
  % max_it          Maximum number of iterations
  % stop_criterion  Stopping criterion for the algorithm
  
  if nargin < 10
    stop_criterion = 1e-10;
  end
  
  for k=1:size(phis,2)
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
      eps_I = 3/16; % inner boundary = 0 because kernel is continuous around 0
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
  scale = 2^(floor(log2((1/4-eps_B)/r_max)));
     

  % Initialize fastsum instances for the different transforms
  kers = cell(1, K);
  inv_kers = cell(1, K);
  for k =1:K
      kers{k} = fastsum(d,kernel,c*scale,flags,n,p,eps_I,eps_B,nn_oversampled,m);
      inv_kers{k} = fastsum(d,kernel,c*scale,flags,n,p,eps_I,eps_B,nn_oversampled,m);
      if k==K
            kers{k}.x =output.*scale; %
            inv_kers{k}.y =output.*scale; %
      else
            kers{k}.x = mus{k+1}.coord.*scale;
            inv_kers{k}.y = mus{k+1}.coord.*scale;
      end
      kers{k}.y = mus{k}.coord.*scale;
      inv_kers{k}.x = mus{k}.coord.*scale;
  end
%tic;

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

  % Iteratively compute maximizing sequence of dual vectors u1, u2, ..., uK
  d = -50;
  d_old = -100;
  iteration = 0;
  %Marginal_error = 100;
  %Marginal_error_old = 50;
  while abs(d-d_old)>stop_criterion  && iteration <max_it %abs(Marginal_error - Marginal_error_old)>stop_criterion
      tic;
      for n_1 = 1:mus{1}.n
        var = phis{2}(:).*betas{1}(:,n_1);
        Marg1 = var'*exp(-sum((mus{2}.coord-mus{1}.coord(n_1,:)).^2,2)/eta);
        phis{1}(n_1) = mus{1}.mu(n_1)./Marg1;
      end
      d_old = d;
      %Marginal_error_old = Marginal_error;
      d = log(phis{1}(:)')*mus{1}.mu(:);
    
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
          d = d + log(phis{k}(:)')*mus{k}.mu(:);
      end

      %Update betas 
      for n_1 = 1:mus{1}.n
          for k=K:-1:2
              if k == K
                  cbeta = exp(-sum((output(n_1,:)-mus{K}.coord).^2,2)/eta);
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

     Marg1_vec = zeros(mus{1}.n,1);
     for n_1 = 1:mus{1}.n
        var = phis{2}(:).*betas{1}(:,n_1);
        Marg1 = var'*exp(-sum((mus{2}.coord-mus{1}.coord(n_1,:)).^2,2)/eta); 
        phis{1}(n_1) = mus{1}.mu(n_1)./Marg1;
        Marg1_vec(n_1) = Marg1;
     end
     d = real(eta*(d-(phis{1}'*Marg1_vec)));
     ds(iteration+1) = real(d);
     fprintf('Iteration %02d,  d = %g\n', iteration, d)

     Marginal_error = 0;
     for k=1:K
         if k==1
             var = (phis{2}.*exp(-pdist2(mus{1}.coord,mus{2}.coord).^2/eta)');
             k_marginal = phis{1}.*reshape(sum(var.*betas{1},1), [mus{1}.n,1]);%.*phis{1};
         else
            k_marginal = phis{k}.*sum(betas{k-1}.*(phis{1}.*gammas{k-1})', 2);
         end
         Marginal_error = Marginal_error + sum(abs(mus{k}.mu-k_marginal));
      end
      %Marginal_errors(iteration+1) = Marginal_error;
      %fprintf('\nIteration %02d,  Error = %g\n', iteration, Marginal_error)
      iteration = iteration+1;
      times(iteration)=toc;
      disp(times(iteration));
      
  end % end while
  averaged_time = mean(times);


end % end function