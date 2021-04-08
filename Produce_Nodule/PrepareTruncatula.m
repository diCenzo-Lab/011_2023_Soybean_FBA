%% Prepare truncatula model 2


TruncaNewNames = importdata('TruncaNamesMapping_of_shared.NewMetNames');
TruncaOldNames = importdata('TruncaNamesMapping_of_shared.OldMetNames');
TruncaNewNames{end+1} = 'HOMO-CIT[c]';
TruncaOldNames{end+1} = 'MNXM722779[c]';
TruncaNewNames{end+1} = 'CPD-3[c]';
TruncaOldNames{end+1} = 'MNXM1026[c]';
TruncaNewNames{end+1} = 'IRON(II)[c]';
TruncaOldNames{end+1} = 'MNXM111[c]';


for i = 1:length(TruncaOldNames)
    
    metindex=findMetIDs(IntegratedTrunca, TruncaOldNames(i));
    
    % _b mets do not formally exist in the mat version of the model. the following
    % if statement avoids this issue
    
    if (metindex > 0)
        IntegratedTrunca.mets(metindex) = TruncaNewNames(i);
    end
    
end

IntegratedTrunca.rxnsformula = printRxnFormula(IntegratedTrunca, IntegratedTrunca.rxns, false);

sol = optimizeCbModel(IntegratedTrunca);

exchangers = ~cellfun(@isempty, regexp(IntegratedTrunca.rxns,'_tx'));
ExportmodelwithOpenBounds = changeRxnBounds(IntegratedTrunca,IntegratedTrunca.rxns(exchangers),10,'u');

ExportmodelwithOpenBounds = IntegratedTrunca;

ExportmodelwithOpenBounds = changeObjective(ExportmodelwithOpenBounds, 'R_BiomassShoot');
sol = optimizeCbModel(ExportmodelwithOpenBounds );

ExportmodelwithOpenBounds.rxnsformula = printRxnFormula(ExportmodelwithOpenBounds , ExportmodelwithOpenBounds.rxns, false);
sol = optimizeCbModel(ExportmodelwithOpenBounds);

