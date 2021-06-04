function CARDAMOM_DEMO
%
%Step 1. Load CBF file
CBF=CARDAMOM_READ_BINARY_FILEFORMAT([getenv('CARDAMOM_DATA_PATH'),'/CARDAMOM_DATA_DRIVERS_EXAMPLE.cbf']);

%Step 2. Run one MCMC chain
%MCO is a structure that contains  MCMC options 
%"niteration" determines how many accepted samples are required
%"printrate" determines how oftern MCMC progress stats are printed (in
%terms of tested samples)
%"samplerate" determines how often the accepted parameter samples are stored
MCO.niterations=10000;
MCO.printrate=1000;
MCO.samplerate=100;
MCO.mcmcid=119;
%CARDAMOM_RUN_MDF runs the MCMC
%All accepted parameters AND all corresponding variables are stored in "CBR"
CBR=CARDAMOM_RUN_MDF(CBF,MCO);
%Step 3. Run CARDAMOM with optimized parameters (This is done internally - at a basic level - inside "CARDAMOM_RUN_MDF")
%CBR=CARDAMOM_RUN_MODEL(CBF,CBR.PARS);

%Step 4. Plot outputs
figure(1);clf
subplot(2,2,1);set(gca,'FontSize',6);
plot(1:72,CBR.GPP);ylabel('GPP [gC/m2/day]');xlabel('Months since Jan 2010');
subplot(2,2,2);
plot(1:72,CBR.FLUXES(:,:,29));ylabel('ET [kgH2O/m2/day]');xlabel('Months since Jan 2010');
subplot(2,2,3);
hist(1-CBR.PARS(:,2));ylabel('Frequency');xlabel('CUE');
subplot(2,2,4);
hist(CBR.PARS(:,17));ylabel('Frequency');xlabel('Leaf C Mass per area [gC/m2]');




disp('CARDAMOM demo has completed')




end
