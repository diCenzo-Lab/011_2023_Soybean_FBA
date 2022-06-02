function NvsNfix
close ALL
changeCobraSolver ('ibm_cplex');

model=readCbModel('Ureide.mat')

[RIPE]=optimizeCbModel(model);
massNON=(RIPE.f).*24/100;

soil=model;

soil = changeRxnBounds(soil,'Root_NH4_tx', 1000000, 'u');
[RIPEY]=optimizeCbModel(soil);
massN=(RIPEY.f).*24/100;
m1=[];N=[];
for n=0:500:10000;
    

model = changeRxnBounds(model,'Root_NH4_tx', n, 'u');
[RIPE]=optimizeCbModel(model);


roger1=find(contains(model.rxns,'Bacteroid_NIT'));

i=0;
if RIPE.f==0
 i=i+1
 n=n+1;
  N=[N,n];
model = changeRxnBounds(model,'Root_NH4_tx', n, 'u');
[RIPE]=optimizeCbModel(model);
roger1=find(contains(model.rxns,'Bacteroid_NIT'));
N1=RIPE.v(roger1);
m1=[m1,N1];

else
N1=RIPE.v(roger1);
m1=[m1,N1];
  N=[N,n];


end
m1(end)

figure(1)
            plot(N/1000,m1*2/1000,'color','b','LineWidth',4);
            hold on, drawnow 
          %  plot(N,massN*ones(size(N)),'--','color','k','LineWidth',4)
            %            hold on, drawnow 
           % plot(N,massNON*ones(size(N)),'--','color','k','LineWidth',4)
            xlabel('Soil ammonium uptake (\mumol/hr/gDW) ','FontSize',50)
            ylabel('N fixation rate (\mumol/hr/gDW)','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
            axis([0 10 0 8])
             set(gcf, 'PaperUnits', 'inches'); 
 x_width=16 ;y_width=15;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('SoilNNfix','-depsc','-loose');



end
