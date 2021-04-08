%% Prepare sino model 2

IntegratedSino = changeObjective(IntegratedSino, 'NIT');
EXlist = IntegratedSino.rxns(strmatch('EX_', IntegratedSino.rxns));
IntegratedSino = changeRxnBounds(IntegratedSino,EXlist,0,'l');

sol = optimizeCbModel(IntegratedSino);
PlantSinoEXlist = importdata('SinoNamesMapping_of_shared.ExListofShared');
PlantSinoEXlist = strrep(PlantSinoEXlist, '_e0', '(e)');
PlantSinoEXlist{end+1} = 'EX_HCIT(e)';
PlantSinoEXlist{end+1} = 'EX_MOBD(e)';
PlantSinoEXlist{end+1} = 'EX_fe2(e)';
PlantSinoEXlist = intersect(PlantSinoEXlist, IntegratedSino.rxns);

sino_test = changeRxnBounds(IntegratedSino, PlantSinoEXlist, -0.1, 'l');
sol = optimizeCbModel(sino_test)

%rename mets according to Metanetx codes (identified throught
%ShareMetsPipeline.sh)

%import new and old names (files written by PrepareMetsSubstitution.pl)
SinoNewNames = importdata('SinoNamesMapping_of_shared.NewMetNames');
SinoOldNames = importdata('SinoNamesMapping_of_shared.OldMetNames');
SinoOldNames{end+1} = 'hcit_e';
SinoNewNames{end+1} = 'MNXM722779__e';
SinoOldNames{end+1} = 'hcit_c';
SinoNewNames{end+1} = 'MNXM722779__c';
SinoOldNames{end+1} = 'mobd_e';
SinoNewNames{end+1} = 'MNXM1026__e';
SinoOldNames{end+1} = 'mobd_c';
SinoNewNames{end+1} = 'MNXM1026__c';
SinoOldNames{end+1} = 'fe2_e';
SinoNewNames{end+1} = 'MNXM111__e';
SinoOldNames{end+1} = 'fe2_c';
SinoNewNames{end+1} = 'MNXM111__c';
SinoOldNames = strrep(SinoOldNames, '_e', '[e]');
SinoNewNames = strrep(SinoNewNames, '__e', '[e]');
SinoNewNames = strrep(SinoNewNames, '__c', '[c]');
sharedMets = importdata('SinoNamesMapping_of_shared.ExListofShared_MNXcode');
sharedMets{end+1} = 'MNXM722779_e0';
sharedMets{end+1} = 'MNXM1026_e0';
sharedMets{end+1} = 'MNXM111_e0';
sharedMets = strrep(sharedMets, '_e0', '[e]');

for i = 1:length(SinoOldNames)
    
    metindex=findMetIDs(IntegratedSino, SinoOldNames(i));
    
    % _b mets do not formally exist in the mat version of the model. the
    % following if statement avoids this issue
    
    if (metindex > 0)
        
        IntegratedSino.mets(metindex) = SinoNewNames(i);
        
    end
    
end

%remove ex reactions (except demand reactions) and then add new EX
%(EXPS) reactions for model cross-talk

AllEx_IntegratedSino = IntegratedSino.rxns(findExcRxns(IntegratedSino));
AllEx_IntegratedSino = AllEx_IntegratedSino(strmatch('EX', AllEx_IntegratedSino));

%Close all EX reactions
IntegratedSino_withEX = changeRxnBounds(IntegratedSino, AllEx_IntegratedSino, 0, 'l');
IntegratedSino_Export = changeRxnBounds(IntegratedSino, AllEx_IntegratedSino, 0, 'l');
%Find EX reactions involving a shared met

AllReactionsWithMNXCompounds = findRxnsFromMets(IntegratedSino_withEX, sharedMets);
AllEXWithMNXCompounds = AllReactionsWithMNXCompounds(~cellfun(@isempty, regexp(AllReactionsWithMNXCompounds,'EX_')));

%Open them and test growth
IntegratedSino_withEX = changeRxnBounds(IntegratedSino_withEX, AllEXWithMNXCompounds, -0.1, 'l');
IntegratedSino_withEX = changeObjective(IntegratedSino_withEX, 'NIT');
sol = optimizeCbModel(IntegratedSino_withEX);

IntegratedSino_Export = changeObjective(IntegratedSino_Export, 'NIT');
sol = optimizeCbModel(IntegratedSino_Export);


%change name to cross-talking EX reactions (from EX_ to EXCT_ and from cpd codes to MNX codes)
noAtpRxns = {'EX_co(e)'; 'EX_co2(e)'; 'EX_h2(e)'; 'EX_h2o(e)'; 'EX_o2(e)'; 'EX_no2(e)'; 'EX_nh4(e)'};
AllEXWithMNXCompoundsNewName = [];
for j = 1:length(AllEXWithMNXCompounds)
    
    formula = printRxnFormula(IntegratedSino_withEX, AllEXWithMNXCompounds{j}, false);
    [mat,tok] = regexp(formula,'^(\w+)\[', 'match', 'tokens');
    mat{:};
    newname = strcat('EXCT_for_',mat{:},'e]');
    newnameB = strcat('EXCT_rev_',mat{:},'e]');
    met1 = strrep(mat{:}, '[', '[e]');
    met2 = strrep(mat{:}, '[', '[c]');
    met_atp = {'ATP[c]'};
    met_adp = {'ADP[c]'};
    met_pi = {'MNXM9[c]'};
    met_proton = {'MNXM1[c]'};
    rxnID = findRxnIDs( IntegratedSino_withEX,AllEXWithMNXCompounds{j});
    IntegratedSino_withEX.rxns{rxnID}  = newname{:};
    AllEXWithMNXCompoundsNewName{end+1} = cell2mat(newname);
    if strmatch(AllEXWithMNXCompounds{j}, noAtpRxns, 'exact')
        if strmatch('EX_nh4(e)', AllEXWithMNXCompounds{j}, 'exact')
            met_proton2 = {'MNXM1[e]'};
            IntegratedSino_Export = addReaction(IntegratedSino_Export,newnameB{:},...
                [met1, met_proton2, met2, met_proton],[-1 -1 1 1],false,0,...
                [],[],[],[],[],[],[],0);
            IntegratedSino_Export = addReaction(IntegratedSino_Export,newname{:},...
                [met2, met_atp, met1, met_adp, met_pi, met_proton],...
                [-1 -0.25 1 0.25 0.25 0.25],false,0,[],[],[],[],[],[],[],0);
        else
            IntegratedSino_Export = addReaction(IntegratedSino_Export,newname{:},...
                [met2, met1],[-1 1],false,0,[],[],[],[],[],[],[],0);
            IntegratedSino_Export = addReaction(IntegratedSino_Export,newnameB{:},...
                [met1, met2],[-1 1],false,0,[],[],[],[],[],[],[],0);
        end
    else
        IntegratedSino_Export = addReaction(IntegratedSino_Export,newname{:},...
            [met2, met_atp, met1, met_adp, met_pi, met_proton],...
            [-1 -0.25 1 0.25 0.25 0.25],false,0,[],[],[],[],[],[],[],0);
        IntegratedSino_Export = addReaction(IntegratedSino_Export,newnameB{:},...
            [met1, met_atp, met2, met_adp, met_pi, met_proton],...
            [-1 -0.25 1 0.25 0.25 0.25],false,0,[],[],[],[],[],[],[],0);
    end
    
end

%% create array for adding exchange reaction to Export model (for debugging purposes)

AllNewMets = [];
for j =  1:length(AllEXWithMNXCompoundsNewName);
    
    formula = printRxnFormula(IntegratedSino_withEX, AllEXWithMNXCompoundsNewName{j}, false);
    [mat,tok] = regexp(formula,'^(\w+)\[', 'match', 'tokens');
    
    met2 = strrep(mat{:}, '[', '[c]');
    AllNewMets{j} = met2{:};
    
end


sol = optimizeCbModel(IntegratedSino_withEX);
sol = optimizeCbModel(IntegratedSino_Export);


CrossTalk_EX_reactions = IntegratedSino_Export.rxns(~cellfun(@isempty, regexp(IntegratedSino_Export.rxns,'EXCT_')));

AllEX_sinoIntegrated =  IntegratedSino_withEX.rxns(findExcRxns(IntegratedSino_withEX));

printRxnFormula(IntegratedSino_withEX, CrossTalk_EX_reactions, false);

sol = optimizeCbModel(IntegratedSino_withEX);


%% close all EX reactions except those that will be open during symbiosis and test growth

TestSinoModel = IntegratedSino_withEX;
AllEx_IntegratedSinoWithEX = TestSinoModel.rxns(findExcRxns(TestSinoModel));
sol = optimizeCbModel(TestSinoModel);
AllEX_reactions = AllEx_IntegratedSinoWithEX(~cellfun(@isempty, regexp(AllEx_IntegratedSinoWithEX,'EX')));

TestSinoModel = changeRxnBounds(TestSinoModel, AllEX_reactions, 0, 'l');
sol = optimizeCbModel(TestSinoModel);

CrossTalk_in_CombinedModel = TestSinoModel.rxns(~cellfun(@isempty, regexp(TestSinoModel.rxns,'EXCT_')));
TestSinoModel = changeRxnBounds(TestSinoModel , CrossTalk_in_CombinedModel, -.1, 'l');
sol = optimizeCbModel(TestSinoModel);

%% prepare model for models integration

IntegratedSino_withEX.rxnsformula = printRxnFormula(IntegratedSino_withEX, TestSinoModel.rxns, false);
IntegratedSino_Export.rxnsformula = printRxnFormula(IntegratedSino_Export, IntegratedSino_Export.rxns,false);
sol = optimizeCbModel(IntegratedSino_Export);


%% testing growth with all the nutrients that will be available inside the plant model
AllLowers = repmat(-.1, 1, length(AllNewMets));
AllUps = repmat(0, 1, length(AllNewMets));
IntegratedSino_ExportDebug = addExchangeRxn(IntegratedSino_Export, AllNewMets, AllLowers, AllUps);
sol = optimizeCbModel(IntegratedSino_ExportDebug);
