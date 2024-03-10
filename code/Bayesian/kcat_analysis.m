
meankcatDists=[];
for i=1:length(kcatDists)
    tic
    meankcatDists=[meankcatDists, mean(kcatDists{i},2)];
end
%meanmeankcatDists=mean(meankcatDists,1);
%medianmeankcatDists=median(meankcatDists,1)

meankcatDists50=[];
for i=1:50:width(meankcatDists)
   meankcatDists50=[meankcatDists50, meankcatDists(:,i)];
end


priorDist = kcatDists{1,1};
priorDist = mean(priorDist');
posteriorDist = kcatDists{1,200};
posteriorDist = mean(posteriorDist');
Prior_vs_Posterior = [priorDist', posteriorDist'];
Foldchange = (posteriorDist./priorDist)';

%% Takes out reactions that are upregulated and has flux. Outputs a table 
%  with reaction index, kcat fold change, reaction name, kcat, flux

RxnsWithFlux = posteriorModel.rxns(find(sol.x(1:2420)));
ecRxnsWithFlux = posteriorModel.ec.rxns(find(ismember(posteriorModel.ec.rxns, RxnsWithFlux)));
UpregulatedRxnsName = posteriorModel.rxnNames(find(ismember(posteriorModel.rxns, ecRxnsWithFlux)));
UpregulatedRxns = posteriorModel.ec.rxns(find(Foldchange>1));
UpregulatedRxnsWithFlux = UpregulatedRxns(find(ismember(UpregulatedRxns, ecRxnsWithFlux)));
UpregulatedRxns_Name = posteriorModel.rxnNames(find(ismember(posteriorModel.rxns, UpregulatedRxnsWithFlux)));
UpregulatedRxns_Foldchange = Foldchange(find(ismember(posteriorModel.ec.rxns, UpregulatedRxnsWithFlux)));
UpregulatedRxns_Kcat = posteriorModel.ec.kcat(find(ismember(posteriorModel.rxns, UpregulatedRxnsWithFlux)));
UpregulatedRxns_Flux = sol.x(find(ismember(posteriorModel.rxns, UpregulatedRxnsWithFlux)));

%ProteinID = posteriorModel.ec.enzymes(find(ismember(posteriorModel.ec.rxns, UpregulatedRxnsWithFlux)));


Output = [UpregulatedRxnsWithFlux, num2cell(UpregulatedRxns_Foldchange), UpregulatedRxns_Name, num2cell(UpregulatedRxns_Kcat), num2cell(UpregulatedRxns_Flux)];


%%
% Multiply kcat with reaction flux to only analyse the kcats of reactions
% that matter.
sol.x=solglc.x+solgly.x+solxyl.x+solhac.x;
PosteriorKcats_with_flux = posteriorDist(find(sol.x(1:2420)))';
PriorKcats_with_flux = priorDist(find(sol.x(1:2420)))';
Foldchange_nonzero = PosteriorKcats_with_flux./PriorKcats_with_flux;
Prior_vs_posterior_nonzero=[PriorKcats_with_flux, PosteriorKcats_with_flux];






%%
% Investigate upregulated kcats, see which reactions they correspond to
UpregulatedKcats = [posteriorModel.ec.rxns(find(Foldchange_nonzero>1)), num2cell(Foldchange_nonzero(find(Foldchange_nonzero>1)))];

%%
figure
violin(log10(Prior_vs_Posterior),{'Prior','Posterior'});
ylabel('log10(kcat)');
xlabel('generation x 50');

figure
violin(log2(Foldchange));
ylabel('log2(Prior vs Posterior kcat fold change)');

%%
figure
grouporder = {'Prior', 'Posterior'};

violin(log10(abs(Prior_vs_Posterior)),'xlabel', grouporder, 'facecolor',[0.03 0.25 0.5], 'edgecolor', [0.03 0.25 0.5], 'facealpha', 0.5);
ylabel('log10(kcat)');
%%
figure
violin(log2(Foldchange_nonzero),'xlabel', grouporder, 'facecolor',[0.03 0.25 0.5], 'edgecolor', 'none', 'facealpha', 0.5);
ylabel('log2(Prior vs Posterior kcat fold change)');
