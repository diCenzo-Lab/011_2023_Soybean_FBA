%% Load data
%addpath /Users/bethanyh/Documents/matlab_toolboxes/cobratoolbox
%addpath /Users/bethanyh/Documents/

addpath /Users/bethanyh/Documents/tigerNEW
addpath(genpath('/Users/bethanyh/Documents/tigerNEW'))

addpath(genpath('/Users/bethanyh/Documents/matlab_toolboxes/cobratoolbox'))
addpath /Users/bethanyh/Documents/Tn-Core-master/Most-Recent-Version
addpath /Users/bethanyh/Documents/FASTCORE_1
changeCobraSolver ('ibm_cplex');

%% adding ATP requirements for transfer reactions

%%
nodulatedPlant=readCbModel('nodulatedPlant.mat');
%nodulatedPlant=readCbModel('WorkingSoy.mat');
biomassrxns
pumpedprotonfix 
ureides
transportATP
nodulatedPlant=tncore_remove_reactions(nodulatedPlant,{'Nodule_R_PROTON_xt','Root_PROTON_tx','Root_PROTON_xt'});

%rna_sino=readcell('bradyexp.xlsx');
rna_sino = table2cell(readtable('bradyexp.txt'));

load('NoduleEXP.mat');
rna_medic=reactionNames;
% rna_medic = table2cell(readtable('Medicago_TPM_values_average.txt')); % assuming thats the mean of expression data for each gene
% rna_sino = table2cell(readtable('Sinorhizobium_TPM_values_average.txt'));
% rna_medic_full = table2cell(readtable('Medicago_TPM_values.txt'));
% rna_sino_full = table2cell(readtable('Sinorhizobium_TPM_values.txt'));
% rna_medic_original = rna_medic;
% rna_sino_original = rna_sino;

gr=nodulatedPlant.grRules(find(contains(nodulatedPlant.rxns,'Nodule_R_CARBODEHYDRAT-RXN_c')));
ru=nodulatedPlant.rules(find(contains(nodulatedPlant.rxns,'Nodule_R_CARBODEHYDRAT-RXN_c')));
RXN0=find(contains(nodulatedPlant.rxns,'Nodule_R_RXN0-5224_c'));
nodulatedPlant.rules(RXN0)=ru;
nodulatedPlant.grRules(RXN0)=gr;
nodulatedPlant = tncore_remove_reactions(nodulatedPlant, {'Nodule_R_CARBODEHYDRAT-RXN_c','Bacteroid_MNXR96123'}); % added MNXR96123 to test if it was needed - remove later 
%% Increase flux
% nodulatedPlant.ub = nodulatedPlant.ub * 1000;
% nodulatedPlant.lb = nodulatedPlant.lb * 1000;
optimizeCbModel(nodulatedPlant)

%% Set model constraints

% Force use of C4-dicarboxylates
sol = optimizeCbModel(nodulatedPlant);
reactions = nodulatedPlant.rxns(strmatch('EXCT_for', nodulatedPlant.rxns));
toRemove = {};
for n = 1:length(reactions)
    pos = findRxnIDs(nodulatedPlant, reactions{n});
    if abs(sol.x(pos)) < 5
        toRemove = vertcat(toRemove, reactions{n});
    end
end
toRemove = setdiff(toRemove, {'EXCT_for_MNXM25[e]';'EXCT_for_MNXM98[e]';'EXCT_for_MNXM93[e]';'EXCT_for_MNXM1[e]'});
nodulatedPlant2 = changeRxnBounds(nodulatedPlant, toRemove, 5, 'u');
nodulatedPlant2 = changeRxnBounds(nodulatedPlant2, {'EXCT_for_MNXM32[e]';'EXCT_for_MNXM29[e]';'EXCT_for_MNXM303[e]';'EXCT_for_MNXM75[e]'}, 5, 'u');
sol2 = optimizeCbModel(nodulatedPlant2)

reactions = nodulatedPlant2.rxns(strmatch('EXCT_for', nodulatedPlant2.rxns));
toRemove = {};
for n = 1:length(reactions)
    pos = findRxnIDs(nodulatedPlant2, reactions{n});
    if abs(sol2.x(pos)) < 5
        toRemove = vertcat(toRemove, reactions{n});
    end
end
toRemove = setdiff(toRemove, {'EXCT_for_MNXM25[e]';'EXCT_for_MNXM98[e]';'EXCT_for_MNXM93[e]';'EXCT_for_MNXM1[e]'});

nodulatedPlant3 = changeRxnBounds(nodulatedPlant2, toRemove, 5, 'u');
sol3 = optimizeCbModel(nodulatedPlant3)
nodulatedPlant = nodulatedPlant3;
nodulatedPlant=changeRxnBounds(nodulatedPlant,'Leave_R_Photon_tx',2920000,'u');
 nodulatedPlant=changeRxnBounds(nodulatedPlant,'EXCT_for_MNXM1026[e]',nodulatedPlant.ub(3518)*2,'u')
 nodulatedPlant=changeRxnBounds(nodulatedPlant,'EXCT_for_MNXM722779[e]',nodulatedPlant.ub(3520)*2,'u');

%% Set growth rate threshold

sol = optimizeCbModel(nodulatedPlant);
growthThresh = 0.99 * sol.f;

%% Prepare the RNA-seq data 

processRNA;

%% Constrain the nodule 

% Turn on tiger
start_tiger('cplex');

% additional code to get cobra_to_tiger to work
reactionsToModify = findRxnsFromGenes(nodulatedPlant, 'Bacteroid_blr2485');
nodulatedPlant.grRules(findRxnIDs(nodulatedPlant, reactionsToModify.gene_Bacteroid_blr2485(:,1))) =  {'( Bacteroid_blr2485 and Bacteroid_blr2486 and ( Bacteroid_bll0282 or Bacteroid_blr3729 ) and ( Bacteroid_bll0283 or Bacteroid_blr3728 ))'};
nodulatedPlant.rules(findRxnIDs(nodulatedPlant, reactionsToModify.gene_Bacteroid_blr2485(:,1))) = {'( x(679) & x(680) & ( x(22) | x(809) ) & ( x(23) | x(808) ))'};
 
% Produce tiger model
nodulatedPlant_tiger = cobra_to_tiger(nodulatedPlant);

% Split the rnaseq files
medicago_rnaseq = cell2mat(medicago_rnaseq_data_model(:,2));
meliloti_rnaseq = cell2mat(meliloti_rnaseq_data_model(:,2));
medicago_genes = medicago_rnaseq_data_model(:,1);
meliloti_genes = meliloti_rnaseq_data_model(:,1);

% Prepare the input for tncore_multiGIMME
fields = {'Medicago'; 'Meliloti'};
expressStruct = struct();
genesStruct = struct();
threshStruct = struct();
expressStruct.Medicago = medicago_rnaseq;
genesStruct.Medicago = medicago_genes;
threshStruct.Medicago = medic_thresh;
expressStruct.Meliloti = meliloti_rnaseq;
genesStruct.Meliloti = meliloti_genes;
threshStruct.Meliloti = sino_thresh;

save('temp_beforeGIMME.mat');

% Run the modified version of gimme
[geneStates, genesOut, sol, tiger, weightStruct] = tncore_multi_gimme(nodulatedPlant_tiger, ...
    fields, expressStruct, genesStruct, threshStruct, 0.99);

save('temp_afterGIMME.mat');

%% Build context specific model in COBRA based on active rxns %% UP TO HERE

buildContextModelNEW;

%% Fix the rules and grRules
% Just in case, may not be necessary

% Change model name
nodulatedPlant = finalNodulatedPlant;

% Fix the grRules string formatting
for n = 1:length(nodulatedPlant.grRules)
    if isstring(nodulatedPlant.grRules{n})
        nodulatedPlant.grRules{n} = str2mat(nodulatedPlant.grRules{n});
    else
        nodulatedPlant.grRules{n} = nodulatedPlant.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(nodulatedPlant.rules)
    if isstring(nodulatedPlant.rules{n})
        nodulatedPlant.rules{n} = str2mat(nodulatedPlant.rules{n});
    else
        nodulatedPlant.rules{n} = nodulatedPlant.rules{n};
    end
end

% Fix the grRules
for n = 1:length(nodulatedPlant.grRules)
    if ~isempty(nodulatedPlant.grRules{n})
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, ' or or', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(nodulatedPlant.grRules{n}, 'or or ', '');
        nodulatedPlant.grRules{n} = strrep(strrep(strrep(nodulatedPlant.grRules{n}, ' or ', 'AAA'), ' or', ''),'AAA', ' or ');
        nodulatedPlant.grRules{n} = strrep(strrep(strrep(nodulatedPlant.grRules{n}, ' or ', 'AAA'), 'or ', ''),'AAA', ' or ');
    end
end

% Fix the rules
for n = 1:length(nodulatedPlant.rules)
    if ~isempty(nodulatedPlant.rules{n})
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, ' | |', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(nodulatedPlant.rules{n}, '| | ', '');
        nodulatedPlant.rules{n} = strrep(strrep(strrep(nodulatedPlant.rules{n}, ' | ', 'AAA'), '| ', ''),'AAA', ' | ');
    end
end

% Fix the grRules string formatting
for n = 1:length(nodulatedPlant.grRules)
    if isstring(nodulatedPlant.grRules{n})
        nodulatedPlant.grRules{n} = str2mat(nodulatedPlant.grRules{n});
    else
        nodulatedPlant.grRules{n} = nodulatedPlant.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(nodulatedPlant.rules)
    if isstring(nodulatedPlant.rules{n})
        nodulatedPlant.rules{n} = str2mat(nodulatedPlant.rules{n});
    else
        nodulatedPlant.rules{n} = nodulatedPlant.rules{n};
    end
end

% Change model name
constrainedNodule = nodulatedPlant;

%% Save and clean workspace

save('allWorkspace.mat');
save('constrainedNodule.mat', 'constrainedNodule');
clear;

