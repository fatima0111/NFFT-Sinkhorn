% Load NFFT
addpath fastsum

addpath Libs
addpath Output_Circle

% fix random number state
rng(4)
output_folder = "./Output_Circle"; % Directory for saving the output data of both algorithms


step  = 4;
K_start = 3;
K_max = 15;
n_times = 1;
etas = [0.1];
max_it = 10;
direct_sum = [false, true];
direct_sum_names = ["NFFT_SINKHORN.....", "SINKHORN......."];
N_t =1;
for eta = etas
    %% laufen
    n_p =700;
    %% Drawing n_p uniform distributed points
    input = sort(unifrnd(-0.5,0.5,n_p,1)); %  reshape(linspace(-0.5,0.5,n_p),[n_p,1]);
    output = input;
    iter = 0;
    k_values = K_start:step:K_max;
    comp_times = zeros(size(k_values));
    comp_times_nfft = zeros(size(k_values));
    for K = k_values
        ds_nfft_sinkhorn = zeros();
        ds_sinkhorn = zeros();
        iter = iter +1;
        %% Sinkhorn algorithm
        %% NFFT_SINKHORN VS SINKHORN
        for j = 1:2
            fprintf(direct_sum_names(j))
            trafo_direct_sum = direct_sum(j);
            inner_time = zeros(N_t,1);%(5);
            phis = cell(1,K);
            mus =cell(1,K);
            %% create discrete measures from images
            for k=1:K
                %% Allocate points
                C= input;
                n = size(C, 1);
                m = zeros(n,1)+ 1/n;
                mu = Measure(m,C);
                mus{k} = mu;
                phi = ones(mu.s)/mu.n;
                phis{k} = phi;
            end
            for t=1:N_t
                tic;
                [us_rec,ds_rec,averaged_times_rec] = NFFT_Sinkhorn_Circle(mus, phis, eta, output, ...
                                                                    trafo_direct_sum, max_it);
                inner_time(t)=toc;
                fprintf('\nTime nfft %02d,  d = %g\n', t, inner_time(t));
            end
            if trafo_direct_sum == false
                  if K==3
                    ds_nfft_sinkhorn = ds_rec;
                  end
                  comp_times_nfft(iter)=averaged_times_rec;
            else
                if K==3
                    ds_sinkhorn= ds_rec;
                end
                comp_times(iter)=averaged_times_rec;
            end
        end
        %% Plot Sinkhorn function for both methods
        if K==3
            iterations_nfft = 1:size(ds_nfft_sinkhorn,2);
            iterations = 1:size(ds_sinkhorn,2);
            %% Save Matrix
            writematrix([iterations_nfft(:) ds_nfft_sinkhorn(:) ds_sinkhorn(:)], ...
            append(output_folder,sprintf("/convergence_rate/Data/convergence_rate_wrt_K_%d_with_n_p_%d_.dat", ...
                                                                                        K, n_p)),'Delimiter',' ');
            figure;
            hold on
            plot(iterations_nfft, ds_nfft_sinkhorn,  '-', color='green');
            plot(iterations, ds_sinkhorn, '-', color='blue');
            legend('NFFT-Sinkhorn', 'Sinkhorn');
            xlabel('Number of iteration');
            ylabel('Sinkhorn function');
            hold off
            prefix = append(output_folder, "/convergence_rate/Plots/", "convergence_rate_wrt_K_");
            name = append(prefix, string(K), "with_n_p_",string(n_p), ".png");
            saveas(gcf,name)
        end
    end 
    %% Plot time complexity with respect to K for both methods
    figure;
    loglog(k_values, comp_times_nfft,  '-', color='green')
    hold on
    loglog(k_values, comp_times, '-', color='blue')
    legend('NFFT-Sinkhorn', 'Sinkhorn')
    xlabel('k');
    ylabel('times');
    hold off
    prefix = append(output_folder, "/complexity/Plots/", "complexity_wrt_K_");
    name = append(prefix, string(K), "_wrt_n_p_", string(n_p), ".png");
    saveas(gcf,name)
    %% Save MAtrix
    writematrix([k_values(:) comp_times_nfft(:) comp_times(:)],  ...
        append(output_folder,sprintf("/complexity/Data/tree_circle_complexity_wrt_K_%d_with_n_p_%d.dat", ...
                                                                            K, n_p)),'Delimiter',' ');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_folder = "./Output_Circle";
    step  = 2500;
    K=3;
    n_p_start = 200;
    n_p_max = 102e2;
    n_times = 1;
    n_p_values = n_p_start:step:n_p_max;
    comp_times = zeros(size(n_p_values));
    comp_times_nfft = zeros(size(n_p_values));
    iter = 0;
    phis = cell(1,K);
    mus =cell(1,K);
    direct_sum = [false, true];
    direct_sum_names = ["NFFT_SINKHORN.....", "SINKHORN......."];
    N_t =1;
    for n_p = n_p_values
        ds_nfft_sinkhorn = zeros();
        ds_sinkhorn = zeros();
        iter = iter +1;
         fprintf("N_P is EQUAL TO%d\n", n_p);
        if n_p < 1000
                max_it = 10;
        elseif n_p >=1000 && n_p < 5000
                max_it = 5;
        else
            max_it = 2;
        end
        for j = 1:2
            fprintf(direct_sum_names(j))
            trafo_direct_sum = direct_sum(j);
            inner_time = zeros(N_t,1);
            %% Drawing uniform samples
            input = sort(unifrnd(-0.5,0.5,n_p,1)); %reshape(linspace(-0.5,0.5,n_p),[n_p,1]);
            output = input;
            for k=1:K
                %% Allocate points
                C= input;
                n = size(C, 1);
                m = zeros(n,1)+ 1/n;
                mu = Measure(m,C);
                mus{k} = mu;
                phi = ones(mu.s)/mu.n;
                phis{k} = phi;
            end
            for t=1:N_t
                tic;
                 [us_rec,ds_rec,averaged_times_rec] = NFFT_Sinkhorn_Circle(mus, phis, eta, output, ...
                                                                            trafo_direct_sum, max_it);
                 inner_time(t)=toc;
                fprintf('\nTime nfft %02d,  d = %g\n', t, inner_time(t));
            end
            if trafo_direct_sum == false
                  comp_times_nfft(iter)=averaged_times_rec;
                  if n_p == 5200
                    ds_nfft_sinkhorn = ds_rec;
                  end
            else
                comp_times(iter)=averaged_times_rec;
                if n_p == 5200
                    ds_sinkhorn= ds_rec;
                end
            end
        end
        if n_p == 5200
            iterations_nfft = 1:size(ds_nfft_sinkhorn,2);
            iterations = 1:size(ds_sinkhorn,2);
            %% save matrix
            writematrix([iterations_nfft(:) ds_nfft_sinkhorn(:) ds_sinkhorn(:)], ...
                append(output_folder,sprintf("/convergence_rate/Data/convergence_rate_wrt_n_p_%d_with_%d.dat" ...
                                                                                       , n_p,K)),'Delimiter',' ');
            %% Plot Sinkhorn function at each iteration
            figure;
            hold on
            plot(iterations_nfft, ds_nfft_sinkhorn,  '-', color='green')
            plot(iterations, ds_sinkhorn, '-', color='blue')
            legend('NFFT-Sinkhorn', 'Sinkhorn')
            xlabel('Number of iteration');
            ylabel('Sinkhorn function');
            prefix = append(output_folder, "/convergence_rate/Plots/", "convergence_rate_wrt_n_p_");
            name = append(prefix, string(n_p), "_with_K_",string(K), ".png");
            saveas(gcf,name)
        end
    end
    %% Plot time complexity with respect to n_p
    figure;
    n_p_values = n_p_start:step:n_p_max;
    loglog(n_p_values, comp_times_nfft,  '-', color='green')
    hold on
    loglog(n_p_values, comp_times, '-', color='blue')
    legend('NFFT-Sinkhorn', 'Sinkhorn')
    xlabel('N');
    ylabel('times');
    hold off
    prefix = append(output_folder, "/complexity/Plots/", "complexity_wrt_n_p_");
    name = append(prefix, string(n_p), "_with_K_",string(K), ".png");
    saveas(gcf,name)
    %% Save MAtrix
    %log_comp_times_nfft = loglog(n_p_values, comp_times_nfft).YData(:);
    %log_comp_times = loglog(n_p_values, comp_times).YData(:);
    writematrix([n_p_values(:) comp_times_nfft(:) comp_times(:)],  ...
        append(output_folder,sprintf("/complexity/Data/1000_circle_complexity_wrt_n_p_%d_with_K_%d.dat", ...
                                                                                   n_p, K)),'Delimiter',' ');
end