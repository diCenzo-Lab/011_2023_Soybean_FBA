function VaryingSR_massbalance
% code to get a window of costs based on a range of S:R ratios (0.09 to
% 0.26) Default value being 0.17 24th Oct 22 
clear 
close ALL
changeCobraSolver ('glpk');

%model=readCbModel('GeorgeUreideMay22.mat')
model=readCbModel('model_with_databases_subs.mat')

 PlantBiomassReactions = find(ismember(model.rxns,{'R_PlantBiomass'}));
 BiomassRootPos = find(ismember(model.mets,'Root_BiomassRoot[c]'));
  BiomassShootPos = find(ismember(model.mets,'Leave_BiomassShoot[c]'));

roop=[];cost=[];perc=[];noo=[];
 for n=0.09:0.01:0.26 
 model.S(BiomassRootPos,PlantBiomassReactions) = -n;
  model.S(BiomassShootPos,PlantBiomassReactions) = -(1-n);
 surfNet(model,model.rxns(1878));
 nin=1-n;
co2adjust=nin-0.83;
co2rxn=find(contains(model.rxns,'Leave_CO2_tx'))
model=changeRxnBounds(model,'Leave_CO2_tx' ,model.ub(co2rxn)*(1+co2adjust),'u');
model.ub(co2rxn)*(1+co2adjust)
1-n
lightrxn=find(contains(model.rxns,'Leave_Photon_tx'))
model=changeRxnBounds(model,'Leave_Photon_tx' ,model.ub(lightrxn)*(1+co2adjust),'u');

% Get maintenance cost reactions
ATPMaintenanceShoot = find(ismember(model.rxns, {'Leave_MNXR96949_c'; 'Leave_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceRoot = find(ismember(model.rxns, {'Root_MNXR105277_c'; 'Root_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceNodule = find(ismember(model.rxns, {'Nodule_MNXR105277_c'; 'Nodule_H+_ATPase_c'; 'Nodule_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceBacteroid = find(ismember(model.rxns, {'Bacteroid_ATPMR'; 'Bacteroid_MNXR153054'}));

% Update maintenance costs
MaintenanceRequirement = 0.6250 / 180 * 32 * 1000 * 1000;
MaintenanceRequirementNodule = MaintenanceRequirement * 0.02;
MaintenanceRequirementBacteroid = 8400 * 0.02 * 0.3 * 1000;
model.lb(ATPMaintenanceShoot) = MaintenanceRequirement * (1-n);
model.lb(ATPMaintenanceRoot) = MaintenanceRequirement * n;
model.lb(ATPMaintenanceNodule) = MaintenanceRequirementNodule * 0.75 * 0.5;
model.lb(ATPMaintenanceBacteroid) = MaintenanceRequirementBacteroid * 0.25 * 0.5;
 ro=optimizeCbModel(model);
 if ro.f~=0
 noo=[noo,n];

 roop=[roop,ro.f*24/1000]
 perc=[perc,abs(((ro.f*24/1000) - 0.0785)/ 0.0785)]
 Nit=find(contains(model.rxns,'Bacteroid_NIT')); %finding reaction number for nitrogenase
N=ro.v(Nit); % finds the flux of nitrogenase reaction
co2=find(contains(model.rxns,'Nodule_CO2_tx')); % finds the reaction number for CO2 flux out of nodule
C=ro.v(co2); % finds the flux of CO2 out of nodule
CN=12.*(abs(C)./abs(N.*2))./14; % calculation for Carbon per nitrogen
cost=[cost,CN];
 else
  n=n-0.01+0.005;   
       noo=[noo,n];
 model.S(BiomassRootPos,PlantBiomassReactions) = -n;
  model.S(BiomassShootPos,PlantBiomassReactions) = -(1-n);
 surfNet(model,model.rxns(1878));


% Get maintenance cost reactions
ATPMaintenanceShoot = find(ismember(model.rxns, {'Leave_MNXR96949_c'; 'Leave_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceRoot = find(ismember(model.rxns, {'Root_MNXR105277_c'; 'Root_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceNodule = find(ismember(model.rxns, {'Nodule_MNXR105277_c'; 'Nodule_H+_ATPase_c'; 'Nodule_ADENOSINETRIPHOSPHATASE-RXN_c'}));
ATPMaintenanceBacteroid = find(ismember(model.rxns, {'Bacteroid_ATPMR'; 'Bacteroid_MNXR153054'}));

% Update maintenance costs
MaintenanceRequirement = 0.6250 / 180 * 32 * 1000 * 1000;
MaintenanceRequirementNodule = MaintenanceRequirement * 0.02;
MaintenanceRequirementBacteroid = 8400 * 0.02 * 0.3 * 1000;
model.lb(ATPMaintenanceShoot) = MaintenanceRequirement * (1-n);
model.lb(ATPMaintenanceRoot) = MaintenanceRequirement * n;
model.lb(ATPMaintenanceNodule) = MaintenanceRequirementNodule * 0.75 * 0.5;
model.lb(ATPMaintenanceBacteroid) = MaintenanceRequirementBacteroid * 0.25 * 0.5;
 ro=optimizeCbModel(model);
 roop=[roop,ro.f*24/1000];
 perc=[perc,abs(((ro.f*24/1000) - 0.0785)/ 0.0785)]
 Nit=find(contains(model.rxns,'Bacteroid_NIT')); %finding reaction number for nitrogenase
N=ro.v(Nit); % finds the flux of nitrogenase reaction
co2=find(contains(model.rxns,'Nodule_CO2_tx')); % finds the reaction number for CO2 flux out of nodule
C=ro.v(co2); % finds the flux of CO2 out of nodule
CN=12.*(abs(C)./abs(N.*2))./14; % calculation for Carbon per nitrogen
cost=[cost,CN];
 end

 end


figure(1)            
            plot(noo,perc*100,'LineWidth',4)
                         hold on, drawnow 
                        plot(noo,((0.302)*ones(size(noo))*100),'-.','color','k','LineWidth',4)
            legend('Varying R:S','R:S=0.17','Location','Best');

            xlabel('Root:Shoot ratio ','FontSize',40)
            ylabel('Relative growth rate cost (%)','FontSize',40)
              set(gca,'LineWidth',2,'FontSize',40)
            %axis([0 60000 0 120])
          % axis([0 10 0 0.08])
             set(gcf, 'PaperUnits', 'inches'); 
 x_width=24 ;y_width=15;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('RGRvaryingRS_withCo2_massbalance','-depsc','-loose');
 


figure(2)            
            plot(noo,cost,'LineWidth',4)
                         hold on, drawnow 
                        plot(noo,((4.1342)*ones(size(noo))),'-.','color','k','LineWidth',4)
            legend('Varying R:S','R:S=0.17','Location','Best');

            xlabel('Root:Shoot ratio ','FontSize',40)
            ylabel('Carbon cost of N fixation (gC/gN)','FontSize',40)
              set(gca,'LineWidth',2,'FontSize',40)
            %axis([0 60000 0 120])
          % axis([0 10 0 0.08])
             set(gcf, 'PaperUnits', 'inches'); 
 x_width=24 ;y_width=15;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('CNvaryingRS_withco2_massbalance','-depsc','-loose');
 
