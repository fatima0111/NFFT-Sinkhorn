addpath Libs
addpath Output_Circle


%% Input
n_p = 400; % Number of support points
input = sort(unifrnd(-0.5,0.5,n_p,1)); % uniform distribution points
%% Simulate SINKHORN and NFFT_SINKHORN for euler flows (circle structure)
K=5;
n_map = 1; % Permutation map between initial position and 
           % final position
d_sums = [false,true];
output_folder = "./Output_Circle";
etas = [0.05]; % Entropy regularization parameter
max_it = 10;   
direct_sum_names = ["NFFT_SINKHORN.....", "SINKHORN......."];
iterations_circle = 1:max_it;
ds_nfft_sinkhorn_circle = zeros(size(etas,2), max_it); % Sinkhorn function values with NFFT-Sinkhorn algorithm
ds_sinkhorn_circle = zeros(size(etas,2), max_it);      % Sinkhorn function values with Sinkhorn algorithm
iteration = 0;
for eta = etas
    iteration = iteration +1;
    ds_nfft_sinkhorn = zeros();
    ds_sinkhorn = zeros();
    us = cell(1,K);
    mus =cell(1,K);
    N_t =1;
    %% create discrete measures from images
    for k=1:K
        %% Allocate points
        C= input;
        if k==1
            output =set_map(n_map, input);
        end
        % Get the number of points
        n = size(C, 1);
        m = zeros(n,1)+ 1/n;
        mu = Measure(m,C);
        mus{k} = mu;
        phi = ones(mu.s)/mu.n;
        phis{k} = phi ;
    end
    for j=1:size(d_sums,2)
        direct_sum = d_sums(j);
        for t=1:N_t
            tic;
            [phis_rec,ds_rec,averaged_times_rec] = NFFT_Sinkhorn_Circle(mus, phis, ...
                                                                eta, output, direct_sum, max_it);
            if direct_sum == false
                ds_nfft_sinkhorn_circle(iteration, :) = ds_rec;
                marg2D_NFFT_Sinkhorn = Compute_Pair_Marginals(mus, phis_rec, eta, output, true);
                %% run test
                Plot_Euler_Marginals(mus, output, eta, marg2D_NFFT_Sinkhorn, true, output_folder);
            else
               ds_sinkhorn_circle(iteration,:)= ds_rec;
               marg2D_Sinkhorn = Compute_Pair_Marginals(mus, phis_rec, eta, output, true);
               Plot_Euler_Marginals(mus, output, eta, marg2D_Sinkhorn, false, output_folder);
            end
            %fprintf('\nTime nfft %02d,  d = %g\n', t, inner_time(t));
        end
    end
    %% Plot Sinkhorn function for both methods
    figure;
    hold on
    plot(iterations_circle, ds_nfft_sinkhorn_circle(iteration,:),  '-', color='green');
    plot(iterations_circle, ds_sinkhorn_circle(iteration,:), '-', color='blue');
    legend('NFFT-Sinkhorn', 'Sinkhorn');
    xlabel('Number of iteration');
    ylabel('Sinkhorn function');
    prefix = append(output_folder, "/convergence_rate/Plots/", "convergence_rate_wrt_K_");
    name = append(prefix, string(K), "with_n_p_",string(n_p), ".png");
    saveas(gcf,name);
end
