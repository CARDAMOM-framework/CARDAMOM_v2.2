% compare the training data & validation data with the CARDAMOM outputs

function CARDAMOM_output_plot(CBR,CBF,vdata)

month=1:192;
figure(1);clf
subplot(3,1,1);set(gca,'FontSize',16);
plotmultilines(1:192,CBR.GPP);
CBF.OBS.GPP(CBF.OBS.GPP==-9999)=NaN;
hold all; plot(month,CBF.OBS.GPP,'ko-','MarkerSize',5);
hold all; plot(vdata(:,1),vdata(:,3),'ro-','MarkerSize',5);
ylabel('GPP [gC/m2/day]');
set(gca,'Xtick',12:24:192,'Xticklabel',2002:2:2016);

subplot(3,1,2);set(gca,'FontSize',16);
plotmultilines(1:192,CBR.NEE);
CBF.OBS.NBE(CBF.OBS.NBE==-9999)=NaN;
hold all; plot(month,CBF.OBS.NBE,'ko-','MarkerSize',5);
hold all; plot(vdata(:,1),vdata(:,4),'ro-','MarkerSize',5);
ylabel('NEE [gC/m2/day]');
set(gca,'Xtick',12:24:192,'Xticklabel',2002:2:2016);

subplot(3,1,3);set(gca,'FontSize',16);
plotmultilines(1:192,CBR.ET);
CBF.OBS.ET(CBF.OBS.ET==-9999)=NaN;
hold all; plot(month,CBF.OBS.ET,'ko-','MarkerSize',5);
hold all; plot(vdata(:,1),vdata(:,5),'ro-','MarkerSize',5);
ylabel('ET [gC/m2/day]');
set(gca,'Xtick',12:24:192,'Xticklabel',2002:2:2016);
end
