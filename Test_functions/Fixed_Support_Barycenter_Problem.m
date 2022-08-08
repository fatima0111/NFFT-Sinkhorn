direct_sum = false; %true; 
eta = 0.005;
max_it = 150;

%% Image interpolation: H-shape tree graph
root = './images/source/';

% Loading dithering images 
C1 = dlmread(append(root,'heart_1000.txt'));
C4 = dlmread(append(root,'tooth_1000.txt'));
C6 = dlmread(append(root,'duck_1000.txt'));
C7 = dlmread(append(root,'redcross_1000.txt'));

% Rescale image measures and remove redundant points
C1 = unique(mod(C1+0.5, 1)-0.5, 'rows');
C4 = unique(mod(C4+0.5, 1)-0.5, 'rows');
C6 = unique(mod(C6+0.5, 1)-0.5, 'rows');
C7 = unique(mod(C7+0.5, 1)-0.5, 'rows');
Cbarycenter = union(union(C1, union(C4,C6, 'rows'), 'rows'), C7, 'rows');

% Set weights of support points to 1/N
m1= ones(size(C1,1),1)/sum(size(C1,1));
m4= ones(size(C4,1),1)/sum(size(C4,1));
m6= ones(size(C6,1),1)/sum(size(C6,1));
m7= ones(size(C7,1),1)/sum(size(C7,1));
mbarycenter = ones(size(Cbarycenter,1),1)/sum(size(Cbarycenter,1));


% Get the number of support points of each marginals
n1 = size(C1, 1);
n4 = size(C4, 1);
n6 = size(C6, 1);
n7 = size(C7, 1);
nbarycenter = size(Cbarycenter,1);

% Marginal weights
lambda =[0.25,0.25, 0.25,0.25]; 

phi1 = ones(n1,1);
phi4 = ones(n4,1);
phi6 = ones(n6,1);
phi7 = ones(n7,1);
phibarycenter = ones(size(Cbarycenter,1),1); %/nbarycenter;

%% create discrete measures from images
keys = [1,2,3,4,5,6,7];
K=size(keys, 2);
constraints = [2,3,5];
phis = {phi1, phibarycenter,phibarycenter, phi4, phibarycenter, phi6, phi7};
ms = {m1(:), mbarycenter(:),mbarycenter(:), m4(:),mbarycenter(:), m6(:), m7(:)};
coords = {C1,  Cbarycenter, Cbarycenter, C4,Cbarycenter, C6, C7};
root = 1;
edges = {[1,2], [2,3], [3,5], [2,4], [5,6], [5,7]};
instance_tree = Tree(1,keys, edges, ms, coords, constraints, lambda);

% direct_sum =true then run Sinkhorn algorithm
% direct_sum = false run NFFT-Sinkhorn algoroithm
[phis_rec,ds_rec,averaged_times_rec] = NFFT_Sinkhorn_Tree(instance_tree, phis, ...
    eta,  direct_sum, max_it);

% Compute barycenters from optimal dual vector phi_k, k=1,4,6,7
% Parameters for fastsum
  d = 2;                % number of dimensions
  kernel = 'gaussian';  % kernel  K(x) = EXP(-x^2/c^2) 
  c = sqrt(eta);        % kernel parameter
  p = 3;                % degree of smoothness of regularization
  flags = 0;            % flags (could be EXACT_NEARFIELD or NEARFIELD_BOXES)
  n = 156;              % expansion degree
  eps_I = 0;            % inner boundary = 0 because kernel is continuous around 0
  if eta > 0.05
    eps_B = 1/16;       % outer boundary
  else
    eps_B = 0;
  end
  m = p;                % cut-off parameter for internal NFFT
  nn_oversampled=2*n;   % oversampling factor for internal NFFT
  
  %% Compute Barycenter
  direct_sum = true;
  % Rescale nodes because fastsum requires them in ball of radius 1/4-eps_B
  r_max = max([instance_tree.nodes(1).marginal.r + instance_tree.nodes(2).marginal.r, ...
      instance_tree.nodes(3).marginal.r+instance_tree.nodes(4).marginal.r]);
  scale = 2^(floor(log2((1/4-eps_B)/r_max)));


  % Initialize fastsum instances for the 4 different transforms
  K21 = fastsum(d,kernel,c*scale/sqrt(lambda(1)),flags,n,p,eps_I,eps_B,nn_oversampled,m);
  K21.x = instance_tree.nodes(1).marginal.coord.*scale; 
  K21.y = instance_tree.nodes(2).marginal.coord.*scale;
  
  K23 = fastsum(d,kernel,c*scale,flags,n,p,eps_I,eps_B,nn_oversampled,m);
  K23.x = instance_tree.nodes(3).marginal.coord.*scale; 
  K23.y = instance_tree.nodes(2).marginal.coord.*scale; 

  K32 = fastsum(d,kernel,c*scale,flags,n,p,eps_I,eps_B,nn_oversampled,m);
  K32.x = instance_tree.nodes(2).marginal.coord.*scale; 
  K32.y = instance_tree.nodes(3).marginal.coord.*scale; 
  
  K24 = fastsum(d,kernel,c*scale/sqrt(lambda(2)),flags,n,p,eps_I,eps_B,nn_oversampled,m);
  K24.x = instance_tree.nodes(4).marginal.coord.*scale;
  K24.y = instance_tree.nodes(2).marginal.coord.*scale;

  K35 = fastsum(d,kernel,c*scale,flags,n,p,eps_I,eps_B,nn_oversampled,m);
  K35.x = instance_tree.nodes(5).marginal.coord.*scale;
  K35.y = instance_tree.nodes(3).marginal.coord.*scale;

  K53 = fastsum(d,kernel,c*scale,flags,n,p,eps_I,eps_B,nn_oversampled,m);
  K53.x = instance_tree.nodes(3).marginal.coord.*scale;
  K53.y = instance_tree.nodes(5).marginal.coord.*scale;

  K56 = fastsum(d,kernel,c*scale/sqrt(lambda(3)),flags,n,p,eps_I,eps_B,nn_oversampled,m);
  K56.x = instance_tree.nodes(6).marginal.coord.*scale;
  K56.y = instance_tree.nodes(5).marginal.coord.*scale;

  K57 = fastsum(d,kernel,c*scale/sqrt(lambda(4)),flags,n,p,eps_I,eps_B,nn_oversampled,m);
  K57.x = instance_tree.nodes(7).marginal.coord.*scale;
  K57.y = instance_tree.nodes(5).marginal.coord.*scale;

  

    K21.alpha = phis_rec{1}(:);
    if direct_sum ==false
        fastsum_trafo(K21); 
    else
        fastsum_trafo_direct(K21);
    end
    a21 = reshape(real(K21.f),instance_tree.nodes(2).marginal.s);
    %a21(a21<0)=0;
    
    K24.alpha = phis_rec{4}(:);
    if direct_sum ==false
        fastsum_trafo(K24); 
    else
        fastsum_trafo_direct(K24);
    end
    a24 = reshape(real(K24.f),instance_tree.nodes(2).marginal.s);
    %a24(a24<0)=0;
    
    K56.alpha = phis_rec{6}(:);
    if direct_sum ==false
        fastsum_trafo(K56); 
    else
        fastsum_trafo_direct(K56);
    end
    a56 = reshape(real(K56.f),instance_tree.nodes(5).marginal.s);
    %a56(a56<0)=0;
    
    
    K57.alpha = phis_rec{7}(:);
    if direct_sum ==false
        fastsum_trafo(K57); 
    else
        fastsum_trafo_direct(K57);
    end
    a57 = reshape(real(K57.f),instance_tree.nodes(5).marginal.s);
    %a36(a36<0)=0;
    
    K35.alpha = a56.*a57;
    if direct_sum ==false
        fastsum_trafo(K35); 
    else
        fastsum_trafo_direct(K35);
    end
    a35 = reshape(real(K35.f),instance_tree.nodes(3).marginal.s);
    %a35(a35<0)=0;
    
    
    K23.alpha = a35;
    if direct_sum ==false
        fastsum_trafo(K23); 
    else
        fastsum_trafo_direct(K23);
    end
    a23 = reshape(real(K23.f),instance_tree.nodes(2).marginal.s);
    
    K32.alpha = a21.*a24;
    if direct_sum ==false
        fastsum_trafo(K32); 
    else
        fastsum_trafo_direct(K32);
    end
    a32 = reshape(real(K32.f),instance_tree.nodes(3).marginal.s);
    
    K53.alpha = a32;
    if direct_sum ==false
        fastsum_trafo(K53);
    else
        fastsum_trafo_direct(K53);
    end
    a53 = reshape(real(K53.f),instance_tree.nodes(5).marginal.s);
    %a53(a53<0)=0;
    
    
    % Compute barycenter mu2
    mbarycenter2 = a21.*a23.*a24;
    mbarycenter2(mbarycenter2<0)=1e-10;
    
    % Compute barycenter mu3
    mbarycenter3 = a32.*a35;
    mbarycenter3(mbarycenter3<0)=1e-10;
    
    % Compute barycenter mu2
    mbarycenter5 = a53.*a56.*a57;
    mbarycenter5(mbarycenter5<0)=1e-10;

    %% Plot figure   
    f2=figure(1);
    subplot(3,3,1)
    scatter(C1(:,1), C1(:,2), m1, '.');
    box on;
    set(gca,'XTick',[], 'YTick', []);
    
    subplot(3,3,2)
    scatter(Cbarycenter(:,1), Cbarycenter(:,2), 1e4*mbarycenter2, '.');
    set(gca,'XTick',[], 'YTick', [])
    box on;
    
    subplot(3,3,3)
    scatter(C4(:,1), C4(:,2), m4, '.')
    set(gca,'XTick',[], 'YTick', []);
    box on;
    
    subplot(3,3,5)
    scatter(Cbarycenter(:,1), Cbarycenter(:,2), 1e4*mbarycenter3, '.');
    set(gca,'XTick',[], 'YTick', [])
    box on;
    
    subplot(3,3,7)
    scatter(C6(:,1), C6(:,2), m6, '.');
    set(gca,'XTick',[], 'YTick', [])
    box on;
      
    subplot(3,3,8)
    scatter(Cbarycenter(:,1), Cbarycenter(:,2), 1e4*mbarycenter5, '.');
    set(gca,'XTick',[], 'YTick', [])
    box on;
    
    subplot(3,3,9)
    scatter(C7(:,1), C7(:,2), m7, '.');
    set(gca,'XTick',[], 'YTick', [])
    box on;