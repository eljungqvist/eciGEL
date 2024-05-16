function variability_comparison(model_conv, ecModel)

[minFluxEc, maxFluxEc] = ecFVA(ecModel, model_conv);
[minFlux, maxFlux] =ecFVA(model_conv,model_conv);

fluxRangeEc = maxFluxEc - minFluxEc;
fluxRange = maxFlux - minFlux;
hold on
cdfplot(fluxRangeEc)
cdfplot(fluxRange)
set(gca, 'XScale', 'log', 'Xlim', [1e-8 1e4])
set(findall(gca, 'Type', 'Line'), 'LineWidth', 2)


end