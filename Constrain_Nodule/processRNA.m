%% Prepare expression data and thresholds

% Remove Medicago genes with no expression
x = 0;
for n = 1:length(rna_medic)
    test = rna_medic{n,2};
    if test ~= 0
        x = x + 1;
        rna_medic_temp(x,:) = rna_medic(n,:);
    end
end
rna_medic = rna_medic_temp;


% Set the expression threshold for Medicago
medic_thresh = mean(mean(cell2mat(rna_medic(:,2:end)))) * 1.1;

% Set the expression threshold for Sinorhizobium
sino_thresh = mean(mean(cell2mat(rna_sino(:,2:end)))) * 1.1;

%% Pre-process the RNA-seq data

% Pull out data for the model plant genes
plantGenes = nodulatedPlant.genes(strmatch('Nodule_', nodulatedPlant.genes));
for n= 1:length(plantGenes)
    plantGenes{n} = strrep(plantGenes{n}, 'Nodule_', '');
end
plant_rnaseq_model = {};
for n = 1:length(plantGenes)
    pos = strmatch(plantGenes{n}, rna_medic(:,1), 'exact');
    if ~isempty(pos)
        if max(cell2mat(rna_medic(pos,2:end))) >= medic_thresh
            plant_rnaseq_model = vertcat(plant_rnaseq_model, rna_medic(pos,:));
        end
    end
end

% Pull out data for the model bacteroid genes
bacterialGenes = nodulatedPlant.genes(strmatch('Bacteroid_', nodulatedPlant.genes));
for n= 1:length(bacterialGenes)
    bacterialGenes{n} = strrep(bacterialGenes{n}, 'Bacteroid_', '');
end
bacterium_rnaseq_model = {};
for n = 1:length(bacterialGenes)
    pos = strmatch(bacterialGenes{n}, rna_sino(:,1), 'exact');
    if ~isempty(pos)
        if max(cell2mat(rna_sino(pos,2:end))) >= sino_thresh
            bacterium_rnaseq_model = vertcat(bacterium_rnaseq_model, rna_sino(pos,:));
        end
    end
end

% Write files to system
plant_rnaseq_model2 = cell2table(plant_rnaseq_model);
bacterium_rnaseq_model2 = cell2table(bacterium_rnaseq_model);
writetable(plant_rnaseq_model2, 'rnaSeqModel_plant.txt', 'Delimiter', '\t', ...
    'WriteVariableNames', false);
writetable(bacterium_rnaseq_model2, 'rnaSeqModel_bacterium.txt', 'Delimiter', '\t', ...
    'WriteVariableNames', false);
%%

% %% Prepare the RNA-seq data
% 
% % Prepare variables for plant in each zone
 medic_zI = horzcat(rna_medic(:,1), rna_medic(:,2));

% 
% % Prepare variables for plant in each zone
 sino_zIId = horzcat(rna_sino(:,1), rna_sino(:,2));


% Rename plant gene names in RNA seq variables
for n = 1:length(medic_zI)
    medic_zI{n,1} = insertBefore(medic_zI{n,1}, 1, 'Nodule_');
end

% Rename bacteria gene names in RNA seq variables
for n = 1:length(sino_zIId)
    sino_zIId{n,1} = insertBefore(sino_zIId{n,1}, 1, 'Bacteroid_');
   
end

% Combine plant RNA-seq variables
medicago_rnaseq_data = medic_zI;

% Combine bacterium RNA-seq variables
meliloti_rnaseq_data = sino_zIId;

%% Force the shoot and root genes to be on

% Identify the shoot and root genes
rootGenes = nodulatedPlant.genes(strmatch('Root_', nodulatedPlant.genes));
shootGenes = nodulatedPlant.genes(strmatch('Leave_', nodulatedPlant.genes));

% Set root and shoot expression above threshold
for n = 1:length(rootGenes)
    rootGenes{n,2} = 2 * medic_thresh;
end
for n = 1:length(shootGenes)
    shootGenes{n,2} = 2 * medic_thresh;
end

% Add root and shoot to the expression list
medicago_rnaseq_data = vertcat(medicago_rnaseq_data, rootGenes);
meliloti_rnaseq_data = vertcat(meliloti_rnaseq_data, shootGenes);
all_rnaseq_data = vertcat(medicago_rnaseq_data, meliloti_rnaseq_data);

%% Identify model RNA-seq data UP TO HERE NOW 

% Pull out data for all model genes. If doesn't exist, set to 0
all_rnaseq_data_model = cell(length(nodulatedPlant.genes), 2);
for n = 1:length(nodulatedPlant.genes)
    pos = strmatch(nodulatedPlant.genes{n}, all_rnaseq_data(:,1), 'exact');
    if isempty(pos)
        all_rnaseq_data_model{n,1} = nodulatedPlant.genes{n,1};
        all_rnaseq_data_model{n,2} = 0;
    else
        all_rnaseq_data_model{n,1} = nodulatedPlant.genes{n,1};
        all_rnaseq_data_model{n,2} = all_rnaseq_data{pos,2};
    end
end

% Split into separate variables for the plant and bacteria
medicago_rnaseq_data_model = {};
meliloti_rnaseq_data_model = {};
for n = 1:length(all_rnaseq_data_model)
    if strmatch('Leave', all_rnaseq_data_model{n,1});
        medicago_rnaseq_data_model = vertcat(medicago_rnaseq_data_model, all_rnaseq_data_model(n,:));
    elseif strmatch('Root', all_rnaseq_data_model{n,1});
        medicago_rnaseq_data_model = vertcat(medicago_rnaseq_data_model, all_rnaseq_data_model(n,:));
    elseif strmatch('Nodule', all_rnaseq_data_model{n,1});
        medicago_rnaseq_data_model = vertcat(medicago_rnaseq_data_model, all_rnaseq_data_model(n,:));
    elseif strmatch('Bacteroid', all_rnaseq_data_model{n,1});
        meliloti_rnaseq_data_model = vertcat(meliloti_rnaseq_data_model, all_rnaseq_data_model(n,:));
    end
end
