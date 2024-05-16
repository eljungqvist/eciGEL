function iGEL_chemostat(ecModel, modelAdapter, qs_max, step)
    
    csource_idx = find(contains(ecModel.rxns, modelAdapter.params.c_source));
    hac_idx = find(contains(ecModel.rxns, 'EX_Acetate'));
    
    ecModel.lb(hac_idx) = 0; %Temp
    
    qs_range = step:step:qs_max;
    
    mu_mod = [];
    qs_mod = [];
    qHac_mod = [];
    
    for i =1:length(qs_range)
        ecModel.lb(csource_idx) = -qs_range(i);
        sol = optimizeCbModel(ecModel);
        mu_mod = [mu_mod; sol.f];
        qs_mod = [qs_mod; -sol.x(csource_idx)];
        qHac_mod = [qHac_mod; sol.x(hac_idx)];
    end

    for i =1:length(qs_range)
        disp(sprintf('%f %f %f',  ...
        qs_mod(i), mu_mod(i), qHac_mod(i)));
    end
    
    %plot(qs_mod, mu_mod)
    
    %plot(qs_mod, qHac_mod)
    plot(mu_mod, qHac_mod)
    
end