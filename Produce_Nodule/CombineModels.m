%  createModel Create a COBRA model from inputs or an empty model
%  structure if no inputs are provided.
%  
%   model = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList,...
%      lowerBoundList,upperBoundList,subSystemList,grRuleList,geneNameList,...
%      systNameList)
%  
%  INPUTS
%   rxnAbrList            List of names of the new reactions
%   rxnNameList           List of names of the new reactions
%   rxnList               List of reactions: format: {'A -> B + 2 C'}
%                         If the compartment of a metabolite is not
%                         specified, it is assumed to be cytoplasmic, i.e. [c]
%                         
%                         
%                         OPTIONAL INPUTS
%   revFlagList           List of reversibility flag (opt, default = 1)
%   lowerBoundList        List of lower bound (Default = 0 or -vMax)
%   upperBoundList        List of upper bound (Default = vMax)
%   subSystemList         List of subsystem (Default = '')
%   grRuleList            List of gene-reaction rule in boolean format (and/or allowed)
%                         (Default = '');
%   geneNameList          List of gene names (used only for translation
%                         from common gene names to systematic gene names)
%   systNameList          List of systematic names
%  
                        
% Make species specific
for n = 1:length(IntegratedSino_Export.mets)
    if strmatch('MNXM', IntegratedSino_Export.mets{n})
        if strfind(IntegratedSino_Export.mets{n}, '[c]')
            IntegratedSino_Export.mets{n} = strcat('Nodule_', IntegratedSino_Export.mets{n});
        else
            IntegratedSino_Export.mets{n} = strcat('Bacteroid_', IntegratedSino_Export.mets{n});
        end
    else
        IntegratedSino_Export.mets{n} = strcat('Bacteroid_', IntegratedSino_Export.mets{n});
    end
end
IntegratedSino_Export.mets(findMetIDs(IntegratedSino_Export, 'Bacteroid_ATP[c]')) = {'Nodule_ATP[c]'};
IntegratedSino_Export.mets(findMetIDs(IntegratedSino_Export, 'Bacteroid_ADP[c]')) = {'Nodule_ADP[c]'};
IntegratedSino_Export.genes = strcat('Bacteroid_', IntegratedSino_Export.genes);
IntegratedSino_Export.rxns = strcat('Bacteroid_', IntegratedSino_Export.rxns);
IntegratedSino_Export.rxnsformula = printRxnFormula(IntegratedSino_Export, IntegratedSino_Export.rxns, false);
IntegratedSino_Export = rmfield(IntegratedSino_Export, 'grRules');
IntegratedSino_Export = tncore_fix(IntegratedSino_Export);
ExportmodelwithOpenBounds.mets = strcat('Nodule_', ExportmodelwithOpenBounds.mets);
ExportmodelwithOpenBounds.mets(findMetIDs(ExportmodelwithOpenBounds, 'Nodule_HOMO-CIT[c]')) = {'Nodule_MNXM722779[c]'};
ExportmodelwithOpenBounds.genes = strcat('Nodule_', ExportmodelwithOpenBounds.genes);
ExportmodelwithOpenBounds.rxns = strcat('Nodule_', ExportmodelwithOpenBounds.rxns);
ExportmodelwithOpenBounds.rxnsformula = printRxnFormula(ExportmodelwithOpenBounds, ExportmodelwithOpenBounds.rxns, false);
ExportmodelwithOpenBounds = tncore_fix(ExportmodelwithOpenBounds);

% Get lists
rxnNameList = vertcat(ExportmodelwithOpenBounds.rxnNames, IntegratedSino_Export.rxnNames);                        
rxnAbrList = vertcat(ExportmodelwithOpenBounds.rxns, IntegratedSino_Export.rxns);
rxnList = vertcat(ExportmodelwithOpenBounds.rxnsformula, IntegratedSino_Export.rxnsformula);
revFlagList = vertcat(ExportmodelwithOpenBounds.rev, IntegratedSino_Export.rev);
lowerBoundList = vertcat(ExportmodelwithOpenBounds.lb, IntegratedSino_Export.lb);
upperBoundList= vertcat(ExportmodelwithOpenBounds.ub, IntegratedSino_Export.ub);
subSystemList= vertcat(ExportmodelwithOpenBounds.subSystems, IntegratedSino_Export.subSystems);
grRuleList = vertcat(ExportmodelwithOpenBounds.grRules, IntegratedSino_Export.grRules);
geneNameList = vertcat(ExportmodelwithOpenBounds.genes, IntegratedSino_Export.genes);

fprintf('\n\n 6. Combining models');
PSmodelOriginal = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList, lowerBoundList, upperBoundList, subSystemList, grRuleList);
PSmodelOriginal.rxns = strrep(PSmodelOriginal.rxns, 'Bacteroid_EXCT_', 'EXCT_' );
PSmodelOriginal = addReaction(PSmodelOriginal, 'Nodule_R_MNXM1026_tx', {'Nodule_MNXM1026[c]'}, 1, 1, -1000, 1000); % Add missing transporter
PSmodelOriginal = addReaction(PSmodelOriginal, 'Nodule_R_MNXM111_tx', {'Nodule_MNXM111[c]'}, 1, 1, -1000, 1000); % Add missing transporter

PSmodelOriginal = changeObjective(PSmodelOriginal, 'Nodule_R_BiomassShoot');
sol = optimizeCbModel(PSmodelOriginal, 'max');

% flux in ex reactions
sol.x(findRxnIDs(PSmodelOriginal, PSmodelOriginal.rxns(findExcRxns(PSmodelOriginal))));

PSmodelOriginal = changeObjective(PSmodelOriginal, 'Bacteroid_NIT');
sol = optimizeCbModel(PSmodelOriginal);
sol.x(findRxnIDs(PSmodelOriginal, PSmodelOriginal.rxns(findExcRxns(PSmodelOriginal))));

PSmodelOriginalMetsChange = PSmodelOriginal;

%% Identify cross-talk reactions between sino and trunca (in original model)

CrossTalk_in_PSmodel = PSmodelOriginal.rxns(~cellfun(@isempty, regexp(PSmodelOriginal.rxns,'EXCT_')));
%printRxnFormula(PSmodelOriginal, CrossTalk_in_PSmodel);


