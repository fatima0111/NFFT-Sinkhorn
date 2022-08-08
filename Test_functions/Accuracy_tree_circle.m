%% Initialisation
n_p = 1000;
K=3;
d_sums = [false, true];
output_folder = "./Output_Tree";
etas = [0.05, 0.1, 0.5]; 
max_it = 10;
direct_sum_names = ["NFFT_SINKHORN.....", "SINKHORN......."];
iterations_circle = 1:max_it;
iterations_tree = 1:max_it;
ds_nfft_sinkhorn_circle = zeros(size(etas,2), max_it);
ds_sinkhorn_circle = zeros(size(etas,2), max_it);
ds_nfft_sinkhorn_tree = zeros(size(etas,2), max_it);
ds_sinkhorn_tree = zeros(size(etas,2), max_it);
iteration = 0;

%% MOTe with tree structured cost function
for eta = etas
    iteration = iteration+1;
    sprintf("Regularization_Parameter_%d\n", eta)
    %% Simulate Sinkkhorn and NFFT-Sinkhorn with respect to the number of marginals K
     fprintf('Number of samples %02d, \n', n_p); 
     instance_tree = Generate_tree(K); %Tree(1,keys, edges, ms, coords);
        
       for j=1:size(d_sums,2)
           disp(direct_sum_names(j));
           direct_sum = d_sums(j);
           [us_rec,ds_rec,averaged_times_rec] = NFFT_Sinkhorn_Tree(instance_tree, us, eta,  direct_sum, max_it);
           if direct_sum == false
               ds_nfft_sinkhorn_tree(iteration,:) = real(ds_rec);
           else
              ds_sinkhorn_tree(iteration,:)= ds_rec;
           end
        end
        
        figure;
        hold on
        plot(iterations_tree, ds_nfft_sinkhorn_tree(iteration,:),  '-', color='green')
        plot(iterations_tree, ds_sinkhorn_tree(iteration,:), '-', color='blue')
        legend('NFFT-Sinkhorn', 'Sinkhorn')
        xlabel('Number of iterations');
        ylabel('Sinkhorn function');
        hold off
        prefix = append(output_folder, "/convergence_rate/Plots/", "convergence_rate_wrt_K_eta");
        name = append(prefix, string(eta), "_with_n_p_",string(n_p), ".png");
        saveas(gcf,name)
end

   %%  MOTe with circle structured cost function
   iteration =0;
   for eta = [0.005] %etas
       iteration = iteration +1;
        ds_nfft_sinkhorn = zeros();
        ds_sinkhorn = zeros();
        us = cell(1,K);
        mus =cell(1,K);
        N_t =1;
        %% create discrete measures from images
        for k=1:K
            %% Allocate points
            input = sort(unifrnd(-0.5,0.5,n_p,1));
            C= input;
            if k==1
                output = input;
            end
            n = size(C, 1);
            m = zeros(n,1)+ 1/n;
            mu = Measure(m,C);
            mus{k} = mu;
            u = ones(mu.s)/mu.n;
            us{k} = u ;
        end
        for j=1:size(d_sums,2)
            direct_sum = d_sums(j);
            for t=1:N_t
                tic;
                [us_rec,ds_rec,averaged_times_rec] = NFFT_Sinkhorn_Circle(mus, us, eta, ...
                                                                        output, direct_sum, max_it);     
                if direct_sum == false
                    ds_nfft_sinkhorn_circle(iteration, :) = ds_rec;
                else
                   ds_sinkhorn_circle(iteration,:)= ds_rec;
                end
            end
        end
    %% Plot Sinkhorn function for fast-sum and methods
    figure;
    hold on
    plot(iterations_circle, ds_nfft_sinkhorn_circle(iteration,:),  '-', color='green')
    plot(iterations_circle, ds_sinkhorn_circle(iteration,:), '-', color='blue')
    legend('NFFT-Sinkhorn', 'Sinkhorn')
    xlabel('Number of iteration');
    ylabel('Sinkhorn function');
    prefix = append(output_folder, "/convergence_rate/Plots/", "convergence_rate_wrt_K_");
    name = append(prefix, string(K), "with_n_p_",string(n_p), ".png");
    saveas(gcf,name)
    end
 %% Save Matrix
 matrix_sinkhorn_tree = zeros(7, size(iterations_tree,2));
 matrix_sinkhorn_tree(1,:)=iterations_tree;
 matrix_sinkhorn_tree(2,:)=ds_nfft_sinkhorn_tree(1,:);
 matrix_sinkhorn_tree(3,:)=ds_nfft_sinkhorn_tree(2,:);
 matrix_sinkhorn_tree(4,:)=ds_nfft_sinkhorn_tree(3,:);
 matrix_sinkhorn_tree(5,:)=ds_sinkhorn_tree(1,:);
 matrix_sinkhorn_tree(6,:)=ds_sinkhorn_tree(2, :);
 matrix_sinkhorn_tree(7,:)=ds_sinkhorn_tree(3,:);



 matrix_sinkhorn_circle = zeros(7, size(iterations_circle,2));
 matrix_sinkhorn_circle(1,:)=iterations_circle;
 matrix_sinkhorn_circle(2,:)=ds_nfft_sinkhorn_circle(1,:);
 matrix_sinkhorn_circle(3,:)=ds_nfft_sinkhorn_circle(2,:);
 matrix_sinkhorn_circle(4,:)=ds_nfft_sinkhorn_circle(3,:);
 matrix_sinkhorn_circle(5,:)=ds_sinkhorn_circle(1,:);
 matrix_sinkhorn_circle(6,:)=ds_sinkhorn_circle(2,:);
 matrix_sinkhorn_circle(7,:)=ds_sinkhorn_circle(3, :);

 writematrix(matrix_sinkhorn_tree, ...
        append("./Output_Tree",sprintf("/convergence_rate/Data/Sinkhorn_iterates_wrt_K_%d_with_n_p_%d_.dat", ...
                                                                                      K, n_p)), 'Delimiter', ' ');
 writematrix(matrix_sinkhorn_circle, ...
        append("./Output_Circle",sprintf("/convergence_rate/Data/Sinkhorn_iterates_wrt_K_%d_with_n_p_%d_.dat", ...
                                                                                        K, n_p)), 'Delimiter', ' ');