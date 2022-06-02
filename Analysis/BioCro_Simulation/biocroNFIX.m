function biocroNFIX
clear 
close ALL


ctrl = readtable('ctrl.csv');
high=readtable('27.09pct_decrease.csv');
seedc=ctrl(:,33);
DVIctrl=ctrl(2077,33);
seedh=high(:,33);
DVIh=high(2077,33);
t=high(:,347);
cc=[];hh=[];
bloo=table2array(seedc);
noo=table2array(seedh);
DVIh=noo(2077);
DVIc=bloo(2077)
for i = 1:length(bloo)      
c=bloo(i)-DVIc;
if c < 0
    c=0;
    cc=[cc,c];
else
cc=[cc,c];
end
h=noo(i)-DVIh;
if h < 0
    h=0;
hh=[hh,h];
else
    hh=[hh,h];
end
end
ctrlseed=cc;
highseed=hh;
finalc=ctrlseed(end)
finalh=highseed(end)
time=table2array(t);
    
figure(3)
            plot(time,ctrlseed,'color','b','LineWidth',5);
            hold on, drawnow 
            plot(time,highseed,'color','k','LineWidth',5)
                        hold on, drawnow 
           % plot(N,massNON*ones(size(N)),'--','color','k','LineWidth',4)
            xlabel('Day of the year (2002) ','FontSize',50)
            ylabel('Yield (Mg/ha)','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
            axis([200 290 0 5])
                       legend('Control','N fix','Location','Best')
             set(gcf, 'PaperUnits', 'inches'); 
 x_width=16 ;y_width=15;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('biocrosimHIGHRGR','-depsc','-loose');
end
