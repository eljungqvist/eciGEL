%% Reactions in central carbon metabolism with strange grRules

% Pyruvate dehydrogenase complex. grRules set or ' or ', generating several
% isozymes when using Gecko. Change to ' and ' and remove weird gene
% associations. 

PDH = 'R00209';
Remove = {'IB49_08600', 'IB49_08610'};
PDH_grRules = {'IB49_08605 and IB49_15580 and IB49_15575 and IB49_03685 and IB49_15570'};
model = changeGrRules(model, PDH, PDH_grRules, true);

Remove_gRules_R00756 = {'IB49_01215'};
model = changeGrRules(model, 'R00756', 'IB49_05585', true);

Remove_gRules_R01518 = {'IB49_17740'};
model = changeGrRules(model, 'R01518', 'IB49_07355 or IB49_03065', true);

Remove_grRules_R00344 = {'IB49_03790'}; %Should only be associated to Acetyl-CoA carboxylase
model = changeGrRules(model, 'R00344', 'IB49_15675', true);

Add_grRules_R00267 = {'IB49_09780'}; % isocitrate dehydrogenase
model = changeGrRules(model, 'R00267', 'IB49_09780', true);

Remove_gRules_R08549 = {'IB49_03685'};
OGDH_grRules = {'IB49_15390 and IB49_15580 and IB49_15395'};
model = changeGrRules(model, 'R08549', OGDH_grRules, true);

Remove_gRules_R01082 = {'IB49_11610'}; % This gene is associated with EC 4.3.1.1 aspartate ammonia lyase instead
model = changeGrRules(model, 'R01082', 'IB49_12670', true);

Remove_grRules_R00742 = {'IB49_18365'};
ACC_grRules = {'IB49_18360 and IB49_03790 and IB49_05590 and IB49_05595 and IB49_03530 and IB49_03795'};
model = changeGrRules(model, 'R00742', ACC_grRules, true);

Acetolactate_synthase = {'IB49_05185 and IB49_05180','IB49_05185 and IB49_05180','IB49_05185 and IB49_05180','IB49_05185 and IB49_05180'};
model = changeGrRules(model, {'R00006', 'R00226', 'R04673', 'R08648'}, Acetolactate_synthase, true);

Pyruvate_kinase={'IB49_05580','IB49_05580','IB49_05580'};
model = changeGrRules(model, {'R00200', 'R01138', 'R01858'}, Pyruvate_kinase, true);

uridylate_kinase={'IB49_16485'};
model = changeGrRules(model, 'R00158', uridylate_kinase, true);

Cellobiose_PTS = {'IB49_13790 and IB49_13800 and IB49_13425'};
model = changeGrRules(model, 'Cellobiose_importABC', Cellobiose_PTS, true);


exportModel(model, '../../../Models/iGEL629', true);
save('../../../Models/iGEL629',"model");

