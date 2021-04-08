%% script performing growth simulations

PSmodelCheck = PSmodelOriginalMetsChange

%set bounds according to biological evidences

%Ammnoium EXCT_MNXM15_e0
%PSmodelCheck = changeRxnBounds(PSmodelCheck, 'EXCT_MNXM15_e0', 0, 'l');
%Succinate EXCT_MNXM25_e0
PSmodelCheck = changeRxnBounds(PSmodelCheck, 'EXCT_for_MNXM25[e]', 0.001326, 'u');
%Malate EXCT_MNXM98_e0
PSmodelCheck = changeRxnBounds(PSmodelCheck, 'EXCT_for_MNXM98[e]', 0.001122, 'u');
%Oxygen EXCT_MNXM4_e0
PSmodelCheck = changeRxnBounds(PSmodelCheck, 'EXCT_for_MNXM4[e]', 1, 'u');
%aspartate
PSmodelCheck = changeRxnBounds(PSmodelCheck, 'EXCT_for_MNXM42[e]', 1, 'u');



TargetReactions = [];
TargetReactions = { 'Bacteroid_NIT';'Nodule_R_BiomassRoot'; 'EXCT_rev_MNXM15[e]'; 'EXCT_for_MNXM25[e]'; ...
    'EXCT_for_MNXM98[e]'; 'EXCT_for_MNXM4[e]'  };
Compound = {'symbiosis'; 'Plant Biomass'; 'Ammonium'; 'Succinate'; 'Malate'; 'Oxygen' };

PSmodelCheck = changeObjective(PSmodelCheck, 'Nodule_R_BiomassRoot');
fprintf('\n\nChecking solution for plant biomass\n\n');
sol = optimizeCbModel(PSmodelCheck, 'max')



PSmodelCheck = changeObjective(PSmodelCheck, 'Bacteroid_NIT');
fprintf('\n\nChecking solution for symbiosis\n\n');
sol = optimizeCbModel(PSmodelCheck)



%% close all cross-talk EX reactions and test symbiosis
PSmodelCheckAllClosed = PSmodelCheck;
PSmodelCheckAllClosed = changeRxnBounds(PSmodelCheckAllClosed, CrossTalk_EX_reactions, 0, 'b');

PSmodelCheckAllClosed = changeObjective(PSmodelCheckAllClosed, 'Bacteroid_NIT');
fprintf('\n\nChecking solution for symbiosis with all cross-talk reactions closed\n\n');
sol = optimizeCbModel(PSmodelCheckAllClosed)


PSmodelCheckAllClosed = changeObjective(PSmodelCheckAllClosed, 'Nodule_R_BiomassRoot');
sol = optimizeCbModel(PSmodelCheckAllClosed)


%Ammnoium EXCT_MNXM15_e0
PSmodelCheckAllClosed = changeRxnBounds(PSmodelCheckAllClosed, 'EXCT_rev_MNXM15[e]', 1000, 'u');
%Succinate EXCT_MNXM25_e0
PSmodelCheckAllClosed = changeRxnBounds(PSmodelCheckAllClosed, 'EXCT_for_MNXM25[e]', 0.001326, 'u');
%Malate EXCT_MNXM98_e0
PSmodelCheckAllClosed = changeRxnBounds(PSmodelCheckAllClosed, 'EXCT_for_MNXM98[e]', 0.001122, 'u');
%Oxygen EXCT_MNXM4_e0
PSmodelCheckAllClosed = changeRxnBounds(PSmodelCheckAllClosed, 'EXCT_for_MNXM4[e]', 1, 'u');
%Aspartate EXCT_MNXM42_e0
PSmodelCheckAllClosed = changeRxnBounds(PSmodelCheckAllClosed, 'EXCT_for_MNXM42[e]', 0.001, 'u');


fprintf('\n\nChecking solution for symbiosis with selected nutrients\n\n');
sol = optimizeCbModel(PSmodelCheckAllClosed) % It is normal for there to be no solution at this step

%table(TargetReactions(:), sol.x(findRxnIDs(PSmodelCheckAllClosed, TargetReactions(:))), Compound(:))
%table(CrossTalk_EX_reactions(:), sol.x(findRxnIDs(PSmodelCheckAllClosed, CrossTalk_EX_reactions(:))))

%% Identify cross-talk reactions between sino and trunca (in model with names changed)

CrossTalk_in_PSmodelMetsChange = PSmodelCheck.rxns(~cellfun(@isempty, regexp(PSmodelCheck.rxns,'EXCT_')));
%printRxnFormula(PSmodelOriginalMetsChange, CrossTalk_in_PSmodelMetsChange);

%% set Objective Functions

PSmodelCheck = changeObjective(PSmodelCheck, 'Bacteroid_NIT');
PSmodelCheck = changeRxnBounds(PSmodelCheck, 'Nodule_R_BiomassRoot', 0.004, 'b');
sol = optimizeCbModel(PSmodelCheck)
   
table(TargetReactions(:), sol.x(findRxnIDs(PSmodelCheck, TargetReactions(:))), Compound(:))

