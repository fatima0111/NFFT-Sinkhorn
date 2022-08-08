function y = set_map(m, input)
    %Permutation map from initial to particles final position
    map1  = @(x) min(2*x+0.5,-2*x+0.5);
    map2  = @(x) mod(x+0.5,0.5);
    map3  = @(x) -x;
    %input2D = dlmread(join(string([root, "initial/initial_", n_p,".txt"]), ''));
    switch m
        case 1
            disp('Permutation map: min(2*x+0.5,0.5-2*x)');
            map = map1;
        case 2
            disp('Permutation map: mod(x+0.5,0.5)');
            map = map2;
        case 3
            disp('Permutation map: -x');
            map = map3;
        otherwise %Default value
            disp('Default permutation map: identity function');
            map = @(x) x;
    end
    y = map(input);

end