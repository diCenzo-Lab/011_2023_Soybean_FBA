function nodule

% Load the models
medicagoModel=readCbModel('soybeangeorge.mat');
% Get truncatula cellular metabolites
trunca_boundary_mets = ~cellfun(@isempty, strfind(medicagoModel.mets,'[c]'));
truncatula_cytosol_mets = medicagoModel.mets(trunca_boundary_mets);
truncatula_cytosol_mets_ok = [];
for i=1:length(truncatula_cytosol_mets)
    NewMetNameTrunca = strrep(truncatula_cytosol_mets(i), '[', '_');
    NewMetNameTrunca = strrep(NewMetNameTrunca, ']', '');
    truncatula_cytosol_mets_ok{end+1,1} = cell2mat(NewMetNameTrunca);
end

% Write truncatula boundary metabolites to file
fileID = fopen('Truncatula_Cytosol.mets','w');
fprintf(fileID,'%s\n' ,truncatula_cytosol_mets_ok{:});
fclose(fileID);

% Extract trunca metabolites from chem_xref
metacyc = table2cell(readtable('chem_xref_metacyc.txt', 'Delimiter', '\t', 'ReadVariableNames', false));
metacyc_trunca= {};
for n = 1:length(medicagoModel.mets);
    met = ['metacyc:' medicagoModel.mets{n}];
    met = strsplit(met, '[');
    met = met{1};
    pos = strmatch(met, metacyc(:,1), 'exact');
    metacyc_trunca = vertcat(metacyc_trunca, metacyc(pos,1:2));
end
for n = 1:length(metacyc_trunca)
    metacyc_trunca{n,1} = strrep(metacyc_trunca{n,1}, 'metacyc:', '');
end
metacyc_trunca2 = {};
for n = 1:length(metacyc_trunca)
    test = [metacyc_trunca{n,1} '['];
    pos = strmatch(test, medicagoModel.mets);
    met = strrep(strrep(medicagoModel.mets{pos(1)}, '[', '_'), ']', '');
    if length(pos) > 1
        for m = 2:length(pos)
            met = [met ';' strrep(strrep(medicagoModel.mets{pos(m)}, '[', '_'), ']', '')];
        end
    end
    metacyc_trunca2{n,1} = metacyc_trunca{n,2};
    metacyc_trunca2{n,2} = met;
end
writecell(metacyc_trunca2, 'MappedCompoundsTrunca', 'Delimiter', '\t');
system('mv MappedCompoundsTrunca.txt MappedCompoundsTrunca');

medicagoModel.rev=[];
for i=1:length(medicagoModel.rxns)
    if medicagoModel.lb(i) < 0
        medicagoModel.rev(i,1)=1;
    else
        medicagoModel.rev(i,1)=0;
    end
end

save('medicagoModel.mat','medicagoModel');

sino=readCbModel('USDA110_model.mat');

% Get sino boundary metabolites (to be used to feed ShareMetsPipeline.sh)
boundary_mets = ~cellfun(@isempty, strfind(sino.mets, '[e]'));
sino_boundary_mets = sino.mets(boundary_mets);
sino_boundary_mets_ok = [];
for i=1:length(sino_boundary_mets)
    NewMetNameSino= strrep(sino_boundary_mets(i), '[', '_');
    NewMetNameSino = strrep(NewMetNameSino, ']', '');
    sino_boundary_mets_ok{end+1,1} = cell2mat(NewMetNameSino);
end

% Write sino boundary metabolites to file
fileID = fopen('Sino_boundary.mets','w'); 
fprintf(fileID,'%s\n' ,sino_boundary_mets{:});
fclose(fileID);

% Extract sino metabolites from chem_xref
bigg = table2cell(readtable('chem_xref_bigg.txt', 'Delimiter', '\t', 'ReadVariableNames', false));
bigg_sino = {};
for n = 1:length(sino.mets);
    met = ['bigg:' sino.mets{n}];
    met = strsplit(met, '[');
    met = met{1};
    pos = strmatch(met, bigg(:,1), 'exact');
    bigg_sino = vertcat(bigg_sino, bigg(pos,1:2));
end
for n = 1:length(bigg_sino)
    bigg_sino{n,1} = strrep(bigg_sino{n,1}, 'bigg:', '');
end
bigg_sino2 = {};
for n = 1:length(bigg_sino)
    test = [bigg_sino{n,1} '['];
    pos = strmatch(test, sino.mets);
    met = strrep(strrep(sino.mets{pos(1)}, '[', '_'), ']', '');
    if length(pos) > 1
        for m = 2:length(pos)
            met = [met ';' strrep(strrep(sino.mets{pos(m)}, '[', '_'), ']', '')];
        end
    end
    bigg_sino2{n,1} = bigg_sino{n,2};
    bigg_sino2{n,2} = met;
end
writecell(bigg_sino2, 'MappedCompoundsSino', 'Delimiter', '\t');
system('mv MappedCompoundsSino.txt MappedCompoundsSino');

sino.rev=[];
for i=1:length(sino.rxns)
    if sino.lb(i) < 0
        sino.rev(i,1)=1;
    else
        sino.rev(i,1)=0;
    end
end
melilotiModel=sino;
save('melilotiModel.mat','melilotiModel');
clear melilotiModel;
%% Rename models
IntegratedTrunca = medicagoModel;
IntegratedSino=sino;

%% Run bash pipeline to map metabolites

BashPipeline

%% Prepare the models

PrepareTruncatula
PrepareSino

%% Combine the models
CombineModels
noduleModel = PSmodelOriginalMetsChange;

finalizeNodule

CheckGrowth
save('allWorkspace.mat');
save('noduleModel.mat','noduleModel');
end