 function NvsMass
clear 
close ALL
changeCobraSolver ('glpk');
%changeCobraSolver ('ibm_cplex');
%model=readCbModel('OGnewestmodelSEPT21.mat');
%model=readCbModel('highRGRmodel.mat')
%model=readCbModel('FinishedHighRGRApril22.mat');
 %model=readCbModel('GeorgeUreideMay22.mat')
model=readCbModel('mass_charge_balanced_model.mat')
 [RIPE]=optimizeCbModel(model);
massNON=(RIPE.f).*24/1000

 

soil=model;rxnFormula1='1 BiomassPlant[c] + 0 BiomassNodule[c] -> BiomassTotal[c]'
soil=addReaction(soil,'R_Biomass',rxnFormula1, [], 0, 0,1000000, 0);
soil=changeObjective(soil,'R_Biomass')
soil = changeRxnBounds(soil,'R_NoduleBiomass', 0, 'b'); 
nod=soil.rxns(find(contains(soil.rxns,'Nodule_')));
bac=soil.rxns(find(contains(soil.rxns,'Bacteroid_')));
nod=vertcat(nod,bac);
for i=1:length(nod)
    soil = changeRxnBounds(soil,nod(i), 0, 'b'); 
end

soil = changeRxnBounds(soil,'Root_NH4_tx', 1000000, 'u');
[RIPEY]=optimizeCbModel(soil);
massN=(RIPEY.f).*24/1000
m1=[];N=[];NN=[];m11=[];
for n=0:500:10000;
%for n=0:50:3600;

model = changeRxnBounds(model,'Root_NH4_tx', n, 'u');
[RIPE]=optimizeCbModel(model);
soil = changeRxnBounds(soil,'Root_NH4_tx', n, 'u');
[ripe1]=optimizeCbModel(soil);
  i=0;

if RIPE.f==0 
  i=i+1
    n=n+1;
  N=[N,n];  NN=[N,n];

model = changeRxnBounds(model,'Root_NH4_tx', n, 'u');
[RIPE]=optimizeCbModel(model);
m=(RIPE.f).*24/1000;
m1=[m1,m];
soil = changeRxnBounds(soil,'Root_NH4_tx', n, 'u');
[ripe1]=optimizeCbModel(soil);
mm=(ripe1.f).*24/1000;
m11=[m11,mm];
elseif  ripe1.f==0
    i=i+1
    n=n+1;
  NN=[NN,n];  N=[N,n];
model = changeRxnBounds(model,'Root_NH4_tx', n, 'u');
[RIPE]=optimizeCbModel(model);
m=(RIPE.f).*24/1000;
m1=[m1,m];
soil = changeRxnBounds(soil,'Root_NH4_tx', n, 'u');
[ripe1]=optimizeCbModel(soil);
mm=(ripe1.f).*24/1000;
m11=[m11,mm];
else
    m=(RIPE.f).*24/1000;
m1=[m1,m];
  N=[N,n];
      mm=(ripe1.f).*24/1000;
m11=[m11,mm];
  NN=[N,n];
end


end
m1(end)
noo=0:10:210;

Ammonium=transpose(N/1000);
UB100=transpose((massN)*ones(size(N)));
LB0=transpose(((massNON)*ones(size(N))));
NfixSoy=transpose(m1);
nonfixSoy=transpose(m11);
maize=transpose(m1+abs(m1-massN)*0.55);
%jeff = cell2table(horzcat(num2cell(Ammonium)))
jeff = cell2table(horzcat(num2cell(Ammonium),num2cell(UB100),num2cell(LB0),num2cell(NfixSoy),num2cell(nonfixSoy),num2cell(maize)));
    writetable(jeff, 'Bacteroid.txt', 'Delimiter', '\t');


figure(1)            
            plot(N/1000,(massN)*ones(size(N)),'--','color','k','LineWidth',4)
             hold on, drawnow 
            plot(N/1000,((massNON)*ones(size(N))),'-.','color','k','LineWidth',4)
             hold on, drawnow
            plot(N/1000,m1,'color','b','LineWidth',4);
          %  hold on, drawnow 

                        plot(N/1000,m11,'color','g','LineWidth',4);
            hold on, drawnow
         %plot(N/1000,(m1+abs(m1-massN)*0.55),'color','r','LineWidth',4);
            hold on, drawnow
               %         legend('100% N','0 % N' ,'Soybean N fixing','Location','Best');

            legend('100% N','0 % N' ,'Soybean N fixing','Soybean non N fixing','Location','Best');
            xlabel('Soil ammonium uptake \mu mol/g/hr ','FontSize',40)
            ylabel(' Relative growth rate g/g/d','FontSize',40)
              set(gca,'LineWidth',2,'FontSize',40)
            %axis([0 60000 0 120])
          % axis([0 10 0 0.08])
             set(gcf, 'PaperUnits', 'inches'); 
 x_width=24 ;y_width=15;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 print('slide','-depsc','-loose');
 
 cost=((massNON-massN)/massN)*100
%print('NewNvmassDec21','-depsc','-loose');
% 
% figure(2)
%             plot(f(noo),noo,'color','b','LineWidth',4);
%             hold on, drawnow 
%             plot(N,massN*ones(size(N)),'--','color','k','LineWidth',4)
%                         hold on, drawnow 
%             plot(N,massNON*ones(size(N)),'--','color','k','LineWidth',4)
%             xlabel('Soil ammonium uptake (\mumol/hr/gDW) ','FontSize',50)
%             ylabel('Soil N (\mu M)','FontSize',50)
%               set(gca,'LineWidth',2,'FontSize',50)
%             %axis([0 6000 0.025 0.057])
%              set(gcf, 'PaperUnits', 'inches'); 
%  x_width=15 ;y_width=6.8;
%  set(gcf, 'PaperPosition', [0 0 x_width y_width]);
% %print('SoilNMASSapril','-depsc','-loose');
%     function ff=f(x)
%         ff=(103.*x)./(216-x);
%     end
end