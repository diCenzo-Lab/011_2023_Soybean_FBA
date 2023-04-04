%% Load the model

model=readCbModel('mass_charge_balanced_model.mat')

h0=find(contains(model.metFormulas,'H0'))
for n=1:length(h0)
    model.metFormulas{h0(n),1}='H2SO3';
    model.metCharges(h0(n),1)=0;
end

%% Generating multiple database IDs for metabolites

% Add new fields
model.metKeggID={};
model.metBiggID={};
model.metChebiID={};
model.metMetacycID={};

% Load the table
metanetxChem = table2cell(readtable('chem_xref.txt', 'Delimiter', '\t','ReadVariableNames', false));

% Extract just the BIGG compound names
biggMets = {};
x = 0;
for n = 1:length(metanetxChem)
    if strmatch('bigg', metanetxChem{n,1})
        x = x + 1;
        biggMets(x,:) = metanetxChem(n,:);
    end
end

% Extract just the Chebi compound names
chebiMets = {};
x = 0;
for n = 1:length(metanetxChem)
    if strmatch('chebi', metanetxChem{n,1})
        x = x + 1;
        chebiMets(x,:) = metanetxChem(n,:);
    end
end

% Extract just the METACYC compound names
metacycMets = {};
x = 0;
for n = 1:length(metanetxChem)
    if strmatch('metacyc', metanetxChem{n,1})
        x = x + 1;
        metacycMets(x,:) = metanetxChem(n,:);
    end
end

% Extract just the KEGG compound names
keggMets = {};
x = 0;
for n = 1:length(metanetxChem)
    if strmatch('kegg', metanetxChem{n,1})
        x = x + 1;
        keggMets(x,:) = metanetxChem(n,:);
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

% Add other identifiers, where possible
for n = 1:length(cpdNames)
    pos = strmatch(cpdNames{n,2}, biggMets(:,2), 'exact');
    if ~isempty(pos)
        model.metBiggID{n,1} = strrep(biggMets{pos(1),1}, 'bigg:','');
    else
        model.metBiggID{n,1} = '';
    end
    pos = strmatch(cpdNames{n,2}, keggMets(:,2), 'exact');
    if ~isempty(pos)
        model.metKeggID{n,1} = strrep(keggMets{pos(1),1}, 'kegg:','');
    else
        model.metKeggID{n,1} = '';
    end
    pos = strmatch(cpdNames{n,2}, metacycMets(:,2), 'exact');
    if ~isempty(pos)
        model.metMetacycID{n,1} = strrep(metacycMets{pos(1),1}, 'metacyc:','');
    else
        model.metMetacycID{n,1} = '';
    end
    pos = strmatch(cpdNames{n,2}, chebiMets(:,2), 'exact');
    if ~isempty(pos)
        model.metChebiID{n,1} = strrep(chebiMets{pos(1),1}, 'chebi:','');
    else
        model.metChebiID{n,1} = '';
    end
end

clearvars -except model

%% Generating multiple database IDs for metabolites

% Add new fields
model.rxnKeggID={};
model.rxnBiggID={};
model.rxnRheaID={};
model.rxnMetacycID={};

% Load the table
metanetxReac = table2cell(readtable('reac_xref.txt', 'Delimiter', '\t','ReadVariableNames', false));

% Extract just the BIGG reaction names
biggReac = {};
x = 0;
for n = 1:length(metanetxReac)
    if strmatch('bigg', metanetxReac{n,1})
        x = x + 1;
        biggReac(x,:) = metanetxReac(n,:);
    end
end

% Extract just the Rhea reaction names
rheaReac= {};
x = 0;
for n = 1:length(metanetxReac)
    if strmatch('rhea', metanetxReac{n,1})
        x = x + 1;
        rheaReac(x,:) = metanetxReac(n,:);
    end
end

% Extract just the METACYC reaction names
metacycReac = {};
x = 0;
for n = 1:length(metanetxReac)
    if strmatch('metacyc', metanetxReac{n,1})
        x = x + 1;
        metacycReac(x,:) = metanetxReac(n,:);
    end
end

% Extract just the KEGG reaction names
keggReac= {};
x = 0;
for n = 1:length(metanetxReac)
    if strmatch('kegg', metanetxReac{n,1})
        x = x + 1;
        keggReac(x,:) = metanetxReac(n,:);
    end
end

% Extract model compound names
rxnNames = cell(length(model.rxns), 3);
for n = 1:length(model.rxns)
    if strmatch('Leave_', model.rxns{n})
        rxnNames{n,1} = 'Leave_';
        rxnNames{n,2} = strrep(model.rxns{n}, 'Leave_', '');
    elseif strmatch('Root_', model.rxns{n})
        rxnNames{n,1} = 'Root_';
        rxnNames{n,2} = strrep(model.rxns{n}, 'Root_', '');
    elseif strmatch('Nodule_', model.rxns{n})
        rxnNames{n,1} = 'Nodule_';
        rxnNames{n,2} = strrep(model.rxns{n}, 'Nodule_', '');
    elseif strmatch('Bacteroid_', model.rxns{n})
        rxnNames{n,1} = 'Bacteroid_';
        rxnNames{n,2} = strrep(model.rxns{n}, 'Bacteroid_', '');
    else
        rxnNames{n,2} = model.rxns{n};
    end
    splitString = strsplit(rxnNames{n,2}, '_');
    if length(splitString) == 2
        rxnNames{n,2} = splitString{1};
        rxnNames{n,3} = splitString{2};
    end
end

for n = 1:length(rxnNames)
    pos = strmatch(rxnNames{n,2}, biggReac(:,2), 'exact');
    if ~isempty(pos)
        model.rxnBiggID{n,1} = strrep(biggReac{pos(1),1}, 'biggR:','');
    else
        model.rxnBiggID{n,1} = '';
    end
    pos = strmatch(rxnNames{n,2}, keggReac(:,2), 'exact');
    if ~isempty(pos)
        model.rxnKeggID{n,1} = strrep(keggReac{pos(1),1}, 'keggR:','');
    else
        model.rxnKeggID{n,1} = '';
    end
    pos = strmatch(rxnNames{n,2}, metacycReac(:,2), 'exact');
    if ~isempty(pos)
        model.rxnMetacycID{n,1} = strrep(metacycReac{pos(1),1}, 'metacycR:','');
    else
        model.rxnMetacycID{n,1} = '';
    end
    pos = strmatch(rxnNames{n,2}, rheaReac(:,2), 'exact');
    if ~isempty(pos)
        model.rxnRheaID{n,1} = strrep(rheaReac{pos(1),1}, 'rheaR:','');
    else
        model.rxnRheaID{n,1} = '';
    end
end

clearvars -except model

%% Save

save('model_with_databases.mat');
