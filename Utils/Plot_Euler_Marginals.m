function Plot_Euler_Marginals(mus, output, eta,  marg2D, isnfft,output_folder)
    % Input               
    % mus                 Instances of class Measure
    % output              Particles position at final time
    % eta                 Entropy regularization parameter
    % marg2D              (1,k)-th Marginal of optimal plan, k=2,...,K 
    % isnfft              Boolean var - True if optimal plan were computed
                          % with NFFT-Sinkhorn and False otherwise
    % output_folder       Parent directory where the plot will be saved
    
    if isnfft == true
        name_suffixe1 = 'nfft_sink_';
    else
        name_suffixe1 = 'sink_';
    end
    name_suffixe2 = append(name_suffixe1, "eta_", string(eta));
    n_p = mus{1}.n;
    K=size(mus,2);
    rows = ceil((K+1)/3); 
    cols = 3; 
    %% Figure
    figure;
    subplot(rows,cols,1);
    input = mus{1}.coord;
    plot(input, input, color="black")
    title('t=0');
    set(gca,'XTick',[], 'YTick', [])
    for k=2:K
        subplot(rows,cols,k);
        pcolor(marg2D{k-1});
        shading flat;
        caxis([1e-7 4e-5]);
        title(join(['t=',string(k-1),"/",string(K)]));
        set(gca,'XTick',[], 'YTick', [])
        colormap gray;
    end
    %colorbar;
    subplot(rows,cols,K+1);
    input = mus{1}.coord;
    plot(input, output, color="black");
    title('t=1');
    set(gca,'XTick',[], 'YTick', [])
    colormap gray
    colormap(flipud(colormap));
    %prefix = append(output_folder, "/pair_marginal/", append("pair_marginal_",name_suffixe2,"n_p_"));
    %name = append(prefix, string(n_p), "_K_", string(K), ".png");
    %saveas(gcf,name)

    %% Save Plot
    figure;
    plot(input, input, color="black");
    set(gca,'XTick',[], 'YTick', []);
    colormap gray;
    colormap(flipud(colormap));
    name = append(output_folder, ...
        sprintf("/euler_plot/%s/%s_t_0_n_p_%d_.png" ...
                                    ,name_suffixe1, name_suffixe2, n_p));
    saveas(gcf,name);
    plot(input, output, color="black");
    shading flat;
    set(gca,'XTick',[], 'YTick', []);
    colormap gray;
    name  = append(output_folder, ...
        sprintf("/euler_plot/%s/%s_t_1_n_p_%d_.png", ...
                                        name_suffixe1, name_suffixe2,  n_p));
    saveas(gcf,name)
    for k=2:K
        pcolor(marg2D{k-1});
        shading flat;
        caxis([1e-7 4e-5]);
        set(gca,'XTick',[], 'YTick', [])
        colormap gray;
        colormap(flipud(colormap));
        name = append(output_folder, ...
        sprintf("/euler_plot/%s/%s_t_%d_n_p_%d_.png", ...
                                    name_suffixe1, name_suffixe2,  k, n_p));
        saveas(gcf,name)
    end

end