%% Global parameters
d_sums = [false, true];
output_folder = "./Output_Tree";

%% Fix n_p and let K varying
etas = [0.1];
max_it = 20;
k_values = 3:1:15;
n_p = 1e4; 
direct_sum_names = ["NFFT_SINKHORN.....", "SINKHORN......."];
for eta = etas
    sprintf("Regularization_Parameter_%d\n", eta)
    %% Simulate Sinkkhorn and Sinkhorn-nfft with respect to the number of marginals K
    comp_times = zeros(size(k_values));
    comp_times_nfft = zeros(size(k_values));
    iter =0;
    for K = k_values
        ds_nfft_sinkhorn = zeros();
        ds_sinkhorn = zeros();
        fprintf('Number of marginals %02d, \n', K);
        iter = iter +1;
        [instance_tree, phis] = Generate_tree(K); 
        
        for j=1:size(d_sums,2)
            disp(direct_sum_names(j));
            direct_sum = d_sums(j);
            [phis_rec,ds_rec,averaged_times_rec] = NFFT_Sinkhorn_Tree(instance_tree, phis, eta,  direct_sum, max_it);
            fprintf('\n');
            if direct_sum == false
                if K == 3
                  ds_nfft_sinkhorn = real(ds_rec);
                end
                comp_times_nfft(1,iter)=averaged_times_rec;
            else
                if K==3
                    ds_sinkhorn= ds_rec;
                end
                comp_times(1, iter)=averaged_times_rec;
            end
        end
        if K ==3
            iterations_nfft_sinkhorn = 1:1:size(ds_nfft_sinkhorn,2);
            iterations_nfft= 1:1:size(ds_sinkhorn, 2);
            %% save second matrix
            writematrix([iterations_nfft(:) ds_nfft_sinkhorn(:) ds_sinkhorn(:)],  ...
                append(output_folder,sprintf("/convergence_rate/Data/convergence_rate_with_K_%d_eta_%g.dat", ...
                                                                                        K,eta)),'Delimiter',' ');
            %% Plot Sinkhorn function for both methods
            figure;

            hold on
            plot(iterations_nfft_sinkhorn, ds_nfft_sinkhorn,  '-', color='green')
            plot(iterations_nfft, ds_sinkhorn, '-', color='blue')
            legend('NFFT-Sinkhorn', 'Sinkhorn')
            xlabel('Number of iterations');
            ylabel('Sinkhorn function');
            hold off
            prefix = append(output_folder, "/convergence_rate/Plots/", "convergence_rate_wrt_K_eta");
            name = append(prefix, string(eta), "_with_n_p_",string(n_p), ".png");
            saveas(gcf,name)
        end
    end
    
    %% Plot Sinkhorn function for both methods
    figure;
    loglog(k_values, comp_times_nfft,  '-', color='green')
    hold on
    loglog(k_values, comp_times, '-', color='blue')
    legend('NFFT-Sinkhorn', 'Sinkhorn')
    xlabel('K');
    ylabel('time (s)');
    hold off
    prefix = append(output_folder, "/complexity/Plots/", "complexity_wrt_K_eta_");
    name = append(prefix, string(eta), "_with_n_p_",string(n_p),".png");
    saveas(gcf,name)
    %% Save MAtrix
    writematrix([k_values(:) comp_times_nfft(:) comp_times(:)], ...
        append(output_folder,sprintf("/complexity/Data/tree_complexity_with_n_p_%d.dat", n_p)),'Delimiter',' ');


    %% Simulate Sinkkhorn and Sinkhorn-nfft with respect to number of support points n_p
    K=10;
    keys = zeros(1,K);
    phis = cell(1, K);
    n_p_start = 200;
    n_p_max = 252e2;
    iter =0;
    step = 5e3;
    n_times = 1;
    n_p_values = n_p_start:step:n_p_max;
    comp_times = zeros(size(n_p_values));
    comp_times_nfft = zeros(size(n_p_values));
    for n_p = n_p_values
        ds_nfft_sinkhorn = zeros();
        ds_sinkhorn = zeros();
        [instance_tree, phis] = Generate_tree(K);
        fprintf('Number of samples %02d, \n', n_p);
        iter = iter +1;
        for j=1:size(d_sums,2)
            sprintf("%s\n",direct_sum_names(j));
            direct_sum = d_sums(j);
            [phis_rec,ds_rec,averaged_times_rec] = NFFT_Sinkhorn_Tree(instance_tree, phis, eta,  direct_sum, max_it);
            fprintf('\n');
            if direct_sum == false
                if n_p == 5200
                  ds_nfft_sinkhorn = real(ds_rec);
                end
                  comp_times_nfft(1,iter)=averaged_times_rec;
                
            else
                if n_p == 5200
                    ds_sinkhorn= ds_rec;
                end
                comp_times(1, iter)=averaged_times_rec;
            end
        end
        if n_p == 5200
            %% Plot Sinkhorn function for both methods
            figure;
            iterations_nfft_sinkhorn = 1:1:size(ds_nfft_sinkhorn,2);
            iterations_nfft= 1:1:size(ds_sinkhorn, 2);
            hold on
            plot(iterations_nfft_sinkhorn, ds_nfft_sinkhorn,  '-', color='green')
            plot(iterations_nfft, ds_sinkhorn, '-', color='blue')
            legend('NFFT-Sinkhorn', 'Sinkhorn')
            xlabel('Number of iterations');
            ylabel('Sinkhorn function');
            hold off
            prefix = append(output_folder, "/convergence_rate/Plots/", "tree_convergence_rate_wrt_n_p_eta_");
            name = append(prefix, string(eta), "with_K_",string(K), ".png");
            saveas(gcf,name)
        end
    end
    
    %% Plot Sinkhorn function for both methods
    figure;
    loglog(n_p_values, comp_times_nfft,  '-', color='green')
    hold on
    loglog(n_p_values, comp_times, '-', color='blue')
    legend('NFFT-Sinkhorn', 'Sinkhorn')
    xlabel('N');
    ylabel('time (s)');
    hold off
    prefix = append(output_folder, "/complexity/Plots/", "tree_complexity_wrt_n_p_eta_");
    name = append(prefix, string(eta), "_with_K_", string(K),".png");
    saveas(gcf,name)
    
    %% Save MAtrix
    log_comp_times_nfft = loglog(n_p_values, comp_times_nfft).YData(:);
    log_comp_times = loglog(n_p_values, comp_times).YData(:);
    writematrix([n_p_values(:) comp_times_nfft(:) comp_times(:)],  ...
        append(output_folder,sprintf("/complexity/Data/tree_complexity_with_K_%d_eta_%g.dat", K,eta)),'Delimiter',' ');
end