%% Load the data

load('constrainedNodule_amide.mat');
model = constrainedNodule;

% Import METANETX reaction database

metanetxRxns = table2cell(readtable('reac_xref.txt','Delimiter', '\t', 'ReadVariableNames', false));

% Import METANETX compound database
metanetxChem = table2cell(readtable('chem_xref.txt', 'Delimiter', '\t','ReadVariableNames', false));

%% Change reaction names

% Extract just the BIGG reaction names
seedRxns = {};
x = 0;
for n = 1:length(metanetxRxns)
    if strmatch('bigg', metanetxRxns{n,1})
        x = x + 1;
        seedRxns(x,:) = metanetxRxns(n,:);
    end
end

% Extract just the METACYC reaction names
metaRxns = {};
x = 0;
for n = 1:length(metanetxRxns)
    if strmatch('metacyc', metanetxRxns{n,1})
        x = x + 1;
        metaRxns(x,:) = metanetxRxns(n,:);
    end
end

% Extract model reaction names
reactionNames = cell(length(model.rxns), 3);
for n = 1:length(model.rxns)
    if strmatch('Leave_', model.rxns{n})
        if strmatch('Leave_R_', model.rxns{n})
        reactionNames{n,1} = 'Leave_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Leave_R_', '');
        else
        reactionNames{n,1} = 'Leave_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Leave_', '');
        end
    elseif strmatch('Root_', model.rxns{n})
        if strmatch('Root_R_', model.rxns{n})
        reactionNames{n,1} = 'Root_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Root_R_', '');
        else 
            reactionNames{n,1} = 'Root_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Root_', '');
        end
    elseif strmatch('Nodule_', model.rxns{n})
        if strmatch('Nodule_R_', model.rxns{n})
        reactionNames{n,1} = 'Nodule_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Nodule_R_', '');
        else
             reactionNames{n,1} = 'Nodule_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Nodule_', '');
        end
            
    elseif strmatch('Bacteroid_', model.rxns{n})
        reactionNames{n,1} = 'Bacteroid_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Bacteroid_', '');
    else
        reactionNames{n,2} = model.rxns{n};
    end
    splitString = strsplit(reactionNames{n,2}, '_');
    if length(splitString) == 2
        reactionNames{n,2} = splitString{1};
        reactionNames{n,3} = splitString{2};
    end
end

% Change reaction names, where possible
for n = 1:length(model.rxns)
    if strmatch('Bacteroid', reactionNames{n,1})
        rxn = ['biggR:' reactionNames{n,2}];
        pos = strmatch(rxn, seedRxns(:,1), 'exact');
        if ~isempty(pos)
            reactionNames{n,2} = seedRxns{pos,2};
        end
    else
        rxn = ['metacycR:' reactionNames{n,2}];
        pos = strmatch(rxn, metaRxns(:,1), 'exact');
        if ~isempty(pos)
            reactionNames{n,2} = metaRxns{pos,2};
        end
    end
end

for n = 1:length(model.rxns)
    rxn = strcat(reactionNames{n,1}, reactionNames{n,2});
    if ~isempty(reactionNames{n,3})
        rxnTemp = strcat(rxn, '_', reactionNames{n,3});
        if strmatch('MNXR', reactionNames{n,2})
            if strmatch(rxnTemp, model.rxns, 'exact')
                model.rxns{n,1} = strcat(rxn, 'b', '_', reactionNames{n,3});
            else
                model.rxns{n,1} = strcat(rxn, '_', reactionNames{n,3});
            end
        else
            model.rxns{n,1} = strcat(rxn, '_', reactionNames{n,3});
        end
    else
        if strmatch('MNXR', reactionNames{n,2})
            if strmatch(rxn, model.rxns, 'exact')
                model.rxns{n,1} = strcat(rxn, 'b');
            else
                model.rxns{n,1} = rxn;
            end
        else
            model.rxns{n,1} = rxn;
        end
    end
end

%% Change metabolite names

% Extract just the BIGG compound names
seedMets = {};
x = 0;
for n = 1:length(metanetxChem)
    if strmatch('bigg', metanetxChem{n,1})
        x = x + 1;
        seedMets(x,:) = metanetxChem(n,:);
    end
end

% Extract just the METACYC compound names
metaMets = {};
x = 0;
for n = 1:length(metanetxChem)
    if strmatch('metacyc', metanetxChem{n,1})
        x = x + 1;
        metaMets(x,:) = metanetxChem(n,:);
    end
end

% Extract model compound names
cpdNames = cell(length(model.mets), 3);
for n = 1:length(model.mets)
    if strmatch('Leave_', model.mets{n})
        cpdNames{n,1} = 'Leave_';
        cpdNames{n,2} = strrep(model.mets{n}, 'Leave_', '');
    elseif strmatch('Root_', model.mets{n})
        cpdNames{n,1} = 'Root_';
        cpdNames{n,2} = strrep(model.mets{n}, 'Root_', '');
    elseif strmatch('Nodule_', model.mets{n})
        cpdNames{n,1} = 'Nodule_';
        cpdNames{n,2} = strrep(model.mets{n}, 'Nodule_', '');
    elseif strmatch('NoduleI_', model.mets{n})
        cpdNames{n,1} = 'NoduleI_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleI_', '');
    elseif strmatch('NoduleIId_', model.mets{n})
        cpdNames{n,1} = 'NoduleIId_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIId_', '');
    elseif strmatch('NoduleIIp_', model.mets{n})
        cpdNames{n,1} = 'NoduleIIp_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIIp_', '');
    elseif strmatch('NoduleIZ_', model.mets{n})
        cpdNames{n,1} = 'NoduleIZ_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIZ_', '');
    elseif strmatch('NoduleIII_', model.mets{n})
        cpdNames{n,1} = 'NoduleIII_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIII_', '');
    elseif strmatch('BacteroidIId_', model.mets{n})
        cpdNames{n,1} = 'BacteroidIId_';
        cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIId_', '');
    elseif strmatch('BacteroidIIp_', model.mets{n})
        cpdNames{n,1} = 'BacteroidIIp_';
        cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIIp_', '');
    elseif strmatch('BacteroidIZ_', model.mets{n})
        cpdNames{n,1} = 'BacteroidIZ_';
        cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIZ_', '');
    elseif strmatch('Bacteroid_', model.mets{n})
        cpdNames{n,1} = 'Bacteroid_';
        cpdNames{n,2} = strrep(model.mets{n}, 'Bacteroid_', '');
    else
        cpdNames{n,2} = model.mets{n};
    end
    splitString = strsplit(cpdNames{n,2}, '[');
    if length(splitString) == 2
        cpdNames{n,2} = splitString{1};
        cpdNames{n,3} = splitString{2};
    end
end

% Change compound names, where possible
for n = 1:length(model.mets)
    if strmatch('Bacteroid', cpdNames{n,1})
        cpd = ['bigg:' cpdNames{n,2}];
        pos = strmatch(cpd, seedMets(:,1), 'exact');
        if ~isempty(pos)
            cpdNames{n,2} = seedMets{pos,2};
        end
    else
        cpd = ['metacyc:' cpdNames{n,2}];
        pos = strmatch(cpd, metaMets(:,1), 'exact');
        if ~isempty(pos)
            cpdNames{n,2} = metaMets{pos,2};
        end
    end
end
for n = 1:length(model.mets)
    cpd = strcat(cpdNames{n,1}, cpdNames{n,2});
    if ~isempty(cpdNames{n,3})
        cpdTemp = strcat(cpd, '_', reactionNames{n,3});
        if strmatch(cpdTemp, model.mets, 'exact')
            model.mets{n,1} = strcat(cpd, 'b', '[', cpdNames{n,3});
        else
            model.mets{n,1} = strcat(cpd, '[', cpdNames{n,3});
        end        
    else
        if strmatch(cpdTemp, model.mets, 'exact')
            model.mets{n,1} = strcat(cpd, 'b');
        else
            model.mets{n,1} = cpd;
        end        
    end
end

%% Remove duplicates

% Get rid of duplicates
[A,B,C] = checkDuplicateRxn(model, 'S');
C = model.rxns(C);
%C=setdiff(C,'Bacteroid_MNXR153054' );
model = tncore_remove_reactions(model, C);
% double checking duplicates
[D,E,F] = checkDuplicateRxn(model, 'S');
F = model.rxns(F);
model = tncore_remove_reactions(model, F);
% Remove any unused genes
model = tncore_remove(model);
optimizeCbModel(model)

% % Rename two metabolites
 %model.mets{118} = 'Bacteroid_MNXM89621b[c]';
 model.mets{975} = 'Nodule_MNXM3224b[c]';
 model.mets{1179} = 'Nodule_MNXM722712b[c]';

% Renaming reactions
% model.rxns{1889}='Root_MNXR109676_c';
% model.rxns{1890}='Leave_MNXR109676_c';
% model.rxns{1898}='Leave_MNXR152255_c';
% model.rxns{1899}='Leave_MNXR144933_c';

%% Ensure maintenance costs are correct (may not be required)

% Get maintenance cost reactions
ATPMaintenanceShoot = find(ismember(model.rxns, {'Leave_MNXR96949_c'; 'Leave_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceRoot = find(ismember(model.rxns, {'Root_MNXR105277_c'; 'Root_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceNodule = find(ismember(model.rxns, {'Nodule_MNXR105277_c'; 'Nodule_H+_ATPase_c'; 'Nodule_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceBacteroid = find(ismember(model.rxns, {'Bacteroid_ATPMR'; 'Bacteroid_MNXR153054'}));

% Update maintenance costs
MaintenanceRequirement = 0.6250 / 180 * 32 * 1000 * 1000;
MaintenanceRequirementNodule = MaintenanceRequirement * 0.02;
MaintenanceRequirementBacteroid = 8400 * 0.02 * 0.3 * 1000;
model.lb(ATPMaintenanceShoot) = MaintenanceRequirement * 0.83;
model.lb(ATPMaintenanceRoot) = MaintenanceRequirement * 0.17;
model.lb(ATPMaintenanceNodule) = MaintenanceRequirementNodule * 0.75 * 0.5;
model.lb(ATPMaintenanceBacteroid) = MaintenanceRequirementBacteroid * 0.25 * 0.5;
optimizeCbModel(model)

%% Finish
finalNodulatedPlant = model;
save('allWorkspace.mat');
save('finalNodulatedPlant_amide.mat', 'finalNodulatedPlant');
clear


