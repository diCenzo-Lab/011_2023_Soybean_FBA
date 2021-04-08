function CombinedModel = BUILDTISSUE()
model=readCbModel('soybeangeorge.mat');
Leafgns=readcell('LEAFEXP1.xlsx');
Rootgns=readcell('ROOTEXP1.xlsx');
model.rev=[];
for i=1:length(model.rxns)
    if model.lb(i) < 0
        model.rev(i,1)=1;
    else
        model.rev(i,1)=0;
    end
end

% Making two models for leaf and root 
Root=model;
Leave=model;

% Remove Unavailable uptakes
% These are light and all Carbon sources for the Root.
disp('Removing Root Reactions')
RootRemovedReacs = find(contains(Root.rxns,{'R_GLC_tx','R_Sucrose_tx','R_Starch_tx','R_Photon_tx'}));
Root = removeRxns(Root,Root.rxns(RootRemovedReacs));

% And any other uptake for the shoot
% Only retain CO2 , O2, and Light (Day)
disp('Removing Importers from Shoot')
LeaveRemovedReacs = find(contains(Leave.rxns,{'R_NH4_tx','R_Ca_tx','R_Fe_tx','R_GLC_tx','R_H2O_tx','R_K_tx','R_Mg_tx','R_NO3_tx','R_PROTON_tx','R_Pi_tx','R_SO4_tx','R_Starch_tx','R_Sucrose_tx','R_TAG_tx'}));
Leave = removeRxns(Leave,Leave.rxns(LeaveRemovedReacs));

% Rename reaction and metabolites to indicate tissue
Root.rxns = strcat('Root_', Root.rxns);
Root.mets = strcat('Root_', Root.mets);
Leave.rxns = strcat('Leave_', Leave.rxns);
Leave.mets = strcat('Leave_', Leave.mets);

disp('Combining Models');
CombinedModel = CombineModelsG(Leave,Root);

ShootUptake = {'GLN','L-ASPARTATE','GLT','NITRATE','L-ALPHA-ALANINE','PROTON','WATER','SULFATE',...
              'MG+2','FE+2','Pi','CO+2','CPD-3'};

RootDelivery = {'ASN','L-ASPARTATE','L-ALPHA-ALANINE','GLN','GLT','SER',...
                'PRO','GLY','THR','VAL','ILE','LEU','LYS','ARG','HIS','WATER',...
                'PROTON','SUCROSE'};
            
disp('Adding Transporters')
NoTransport = CombinedModel.rxns;
CombinedModel = addReaction(CombinedModel,'TRS_AMMONIUM',{'Root_AMMONIA[c]','Root_PROTON[c]','Leave_PROTON[c]','Leave_AMMONIA[c]'},[-1 -1 1 1],0,0,1000);
for i=1:numel(ShootUptake)
    CombinedModel = addReaction(CombinedModel,['TRS_' ShootUptake{i}],{['Root_' ShootUptake{i} '[c]'],['Leave_' ShootUptake{i} '[c]'] }, [-1 1], 0,0,1000);
end
for i=1:numel(RootDelivery)
    CombinedModel = addReaction(CombinedModel,['TSR_' RootDelivery{i}],{['Leave_' RootDelivery{i} '[c]'],['Root_' RootDelivery{i} '[c]'] }, [-1 1], 0,0,1000);
end
Transporters = setdiff(CombinedModel.rxns,NoTransport);
% Determine 2 Models
% Day Model: Allow Light, No Starch, BiomassShoot for Leave, BiomassRoot for Root
% For Day use Ammonia and Nitrate for growth.          
LightImport = find(contains(CombinedModel.rxns,'Leave_R_Photon_tx'));
RootBiomass = find(contains(CombinedModel.rxns,'Root_R_BiomassRoot'));
LeaveBiomass = find(contains(CombinedModel.rxns,'Leave_R_BiomassShoot'));
NitrateUptake = find(contains(CombinedModel.rxns,'Root_R_NO3_tx'));
Ammoniumuptake = find(contains(CombinedModel.rxns,'Root_R_NH4_tx'));

% loop to generate table of expression data switched on & off against
% GeneIDs 
LEAFY=Leafgns;
for i=1:length(Leafgns)
    if Leafgns{i,2}<21
       LEAFY{i,2}=0;
       LEAFY{i,1}=Leafgns(i,1);
    elseif Leafgns{i,2}==21 || Leafgns{i,2}>21
        LEAFY{i,2}=1;
         LEAFY{i,1}=Leafgns{i,1};
    end
end
GeneID=LEAFY(:,1);
bob=LEAFY(:,2);
bob1=string(bob);
bob2=double(bob1);
[A,B] = ismember(string(GeneID),model.genes);

%Determine the Leave Present Reactions
LeaveActivities = zeros(1,numel(model.genes));
LeaveActivities(B(B~=0)) = bob2(A);
LeaveOn = model.genes(LeaveActivities==1);
LeaveReacs = findRxnsActiveWithGenes(model,LeaveOn);
LeaveReacs = strcat('Leave_',LeaveReacs);

%Determine the Root Present Reactions
ROOTY=Rootgns;
for j=1:length(Rootgns)
    if Rootgns{j,2}<21
       ROOTY{j,2}=0;
       ROOTY{j,1}=Rootgns(j,1);
    elseif Rootgns{j,2}==21 || Rootgns{j,2}>21
        ROOTY{j,2}=1;
         ROOTY{j,1}=Rootgns{j,1};
    end
end
GeneID=ROOTY(:,1);
seb=ROOTY(:,2);
seb1=string(seb);
seb2=double(seb1);
[C,D] = ismember(string(GeneID),model.genes);
RootActivities(D(D~=0)) = seb2(C);
RootOn = model.genes(RootActivities==1);
RootReacs = findRxnsActiveWithGenes(model,RootOn);
RootReacs = strcat('Root_',RootReacs);

% Determine core reactions
CoreReacs = union(RootReacs,LeaveReacs);
DayAmmoniaReacs = union(CoreReacs,CombinedModel.rxns([Ammoniumuptake,LightImport,LeaveBiomass,RootBiomass]));
DayNitrateReacs = union(CoreReacs,CombinedModel.rxns([NitrateUptake,LightImport,LeaveBiomass,RootBiomass]));

% Make consistent model
disp('Making Model Consistent')
ConsistentCombined = fastcc(CombinedModel,1e-4);
ConsistentModel = removeRxns(CombinedModel,setdiff(CombinedModel.rxns,CombinedModel.rxns(ConsistentCombined)));

% Make core models
disp('Generating Models for different Conditions')

StarchImport = find(ismember(ConsistentModel.rxns,'Leave_R_Starch_tx'));
NitrateUptake = find(ismember(ConsistentModel.rxns,'Root_R_NO3_tx')); 
Ammoniumuptake = find(ismember(ConsistentModel.rxns,'Root_R_NH4_tx'));

DayAmmoniumModel = ConsistentModel;
DayAmmoniumModel.ub([StarchImport,NitrateUptake]) = 0;

DayNitrateModel = ConsistentModel;
DayNitrateModel.ub([StarchImport,Ammoniumuptake]) = 0;

disp('Making Models Consistent')
AA = fastcc(DayAmmoniumModel,1e-4);
DayAmmoniumModel = removeRxns(DayAmmoniumModel,setdiff(DayAmmoniumModel.rxns,DayAmmoniumModel.rxns(AA)));
AA = fastcc(DayNitrateModel,1e-4);
DayNitrateModel = removeRxns(DayNitrateModel,setdiff(DayNitrateModel.rxns,DayNitrateModel.rxns(AA)));

CoreDayAmmonium = find(ismember(DayAmmoniumModel.rxns,DayAmmoniaReacs));
NoPenaltyDayAmmonium = find(ismember(DayAmmoniumModel.rxns,Transporters));
CoreDayNitrate = find(ismember(DayNitrateModel.rxns,DayNitrateReacs));
NoPenaltyDayNitrate = find(ismember(DayNitrateModel.rxns,Transporters));

disp('Using fastcore to restricte Day Ammonium Model')
AA = fastcore(CoreDayAmmonium, DayAmmoniumModel,1e-4);
DayAmmoniumModel = removeRxns(DayAmmoniumModel,setdiff(DayAmmoniumModel.rxns,DayAmmoniumModel.rxns(AA)));

disp('Using fastcore to restricte Day Nitrate Model')
AA = fastcore(CoreDayNitrate,DayNitrateModel,1e-4);
DayNitrateModel = removeRxns(DayNitrateModel,setdiff(DayNitrateModel.rxns,DayNitrateModel.rxns(AA)));

disp('Adding Combined Biomass Reaction')

%Combine All Models
CombinedModelReacs = union(Transporters,union(DayAmmoniumModel.rxns,DayNitrateModel.rxns));

%Keep the Root Proton Exchanger
CombinedModelReacs{end+1} = 'Root_R_PROTON_tx';
CombinedModel = removeRxns(CombinedModel,setdiff(CombinedModel.rxns,CombinedModelReacs));

%Now, add a reaction combining Biomass for root and shoot
CombinedModel = addReaction(CombinedModel,'R_Biomass',{'Leave_BiomassShoot[c]','Root_BiomassRoot[c]','Biomass[c]'} ,[-1 -1 1],0,0,1000);
CombinedModel = addReaction(CombinedModel,'EX_Biomass',{'Biomass[c]'} ,[-1],0,0,1000);
ShootBiomassReactions = find(ismember(CombinedModel.rxns,{'Leave_R_BiomassShoot'}));
RootBiomassReactions = find(ismember(CombinedModel.rxns,{'Root_R_BiomassRoot'}));
BiomassShootPos = find(ismember(CombinedModel.mets,'Leave_BiomassShoot[c]'));
BiomassRootPos = find(ismember(CombinedModel.mets,'Root_BiomassRoot[c]'));
CombinedModel.S(BiomassShootPos,ShootBiomassReactions) = 1;
CombinedModel.S(BiomassRootPos,RootBiomassReactions) = 1;
CombinedModel.S(BiomassShootPos,end-1) = -2/3;
CombinedModel.S(BiomassRootPos,end-1) = -1/3;
CombinedModel.c(:) = 0;
CombinedModel.c(end) = 1;
CombinedModel.lb(findRxnIDs(CombinedModel, 'Root_EX_BiomassRoot')) = 0;
CombinedModel.lb(findRxnIDs(CombinedModel, 'Leave_EX_BiomassShoot')) = 0;

end