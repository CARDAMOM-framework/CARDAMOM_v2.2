function MCMCoutput=MHMCMC_EXAMPLE
%*************Example code Summary*************
%
%The following example code demonstrates how the MHMCMC.m matlab code can be used
%to combine a user-defined model and model observational constraints to estimate model parameter probability distribution.
%For more details on the Metropolis-Hastings Markov Chain Monte Carlo (MHMCMC) approach and implementation, see
%Bloom & Williams, 2015; Ziehn et al., 2012; Xu et al., 2006, and
%references therein; references are at the end of this code).
%
%Specifically, this example shows how a number of model parameters (p1,p2,p3,...) can be estimated
%based on a set of observations (o1, o2, o3,...) and their associated uncertainties 
%(u1, u2, u3,...). All necessary information (model, data, parameter info)
%are externally defined (as shown in this example code) and are passed to
%MHMCMC.m. The MHMCMC.m routine returns anumber (e.g. 10,000) parameter
%"samples"; these can then be used to derive the parameter probability
%distribution.
%
%The following model (a wetland CH4 emissions model which varies
%with temperature) will be used as an example.
%
%MCH4 = p1 * p2.^(T/10),
%
%where MCH4 are monthly model CH4 emissions, and T are monthly mean temperature values
%between january and december. The unknown model parameter values p1 and p2
%will be derived by minimizing the difference between model CH4 emissions
%MCH4 to monthly observations OCH4 (defined in STEP 1).
%
%To use the MHMCMC.m code for estimation of parameters p1 and p2, the
%following inputs are required: 
%(1) DATA: a structure containing all model drivers and observations
%(2) MLF: the model-likelihood function handle
%(3) PARS: a structure defining the prior information on parameters p1 and p2
%(4) MCO: a structure containing MHMCMC execution options, such as number of
%samples required (N), etc.
%
%The MHMCMC inputs are described in more detail within the example code.
%The exectution of the MHMCMC code will return an output structure,
%MCMCoutput:
%
%MCMCoutput=MHMCMC(DATA,MLF,PARS,MCO);
%
%This structure contains N parameter samples, as well as additional information on
%the MHMCMC run. Based on the N parameter samples, a probability density
%function (e.g. a histogram) of each parameter can be derived.
%
%Note: This is example code can be used as a template for a new project
%involving a MHMCMC optimization approach. The MHMCMC.m code has been
%extensively tested, and has been written such that no changes to MHMCMC.m are
%necessary for each new application (e.g. new model, new observations, etc.).
%For most MHMCMC applications, all changes can simply be made in the MHMCMC_EXAMPLE code. 
%
%Function last edited on Aug 6th 2015 by A. A. Bloom
%abloom@jpl.nasa.gov







%*************STEP 1. Store all datasets in DATA structure. *************
%Example model drivers and observations:
%
%Monthly Temperature  (model driver)
T=[6,5,7,10,12,15,19,22,20,17,11,6];
%Monthly observations: CH4 fluxes over wetland area (mg m-2 day-1)
OCH4=[ 13,16,15,22,21,29,39,45,42,31,19,16];
%Monthly Observation uncertainty (mg m-2 day-1).
UCH4=[1,2,1,2,1,2,1,2,1,2,1,2];
%All data stored here:
DATA.temp=T;
DATA.obs=OCH4;
DATA.obsunc=UCH4;

%*************STEP 2. Define model likelihood function MLF*************
%INPUTS: 
%   DATA: all data (DATA structure)
%   pars: vector with parameter values (in example case, pars(1) is
%   equivalent to p1, and pars(2) to p(2)).
%OUTPUTS: 
%   P = log probability (range = -inf to 0)
%
%MLF must be defined as a function handle (see below). For example purposes, the full model and model-likelihood functions have already
%been defined (see MODEL_LIKELIHOOD_FUNCTION and WETLAND_CH4_MODEL functions at the bottom of this script)
%
%Example:
MLF=@(DATA,pars) MODEL_LIKELIHOOD_FUNCTION(DATA,pars);

%Test Step:
%
%To test MLF, try running with mock parameter values
%       pars_example=[10,2];
%       P=MLF(DATA,pars_example);
%To test the model only, try running WETLAND_CH4_MODEL
%       pars_example=[10,2];
%       MCH4 = WETLAND_CH4_MODEL(DATA,pars);
%
%To try the above tests before the MHMCMC (not necessary), un-comment out the next three lines 
% disp('Test Step: model and MLF functions (see Step 2 instructions)')
% disp('When done, type "return"')
% keyboard

%*************STEP 3. Define PARS structure*************
%
%The PARS structure contains all parameter information needed to find the probability
%distribution of parameter values p1 and p2.
%The required fields are the minimum and maximum values of p1 and p2
%Optional fields include the initial MCMC step size (default = 0.001) and the initial values
%of p1 and p2 (default = random). 
%
%PARS fields:
%
%PARS.min = minimum values
%PARS.max = maximum values
%PARS.step = initial step size (optional)
%PARS.init = initial values (optional)
%
%Example:
PARS.min=[0.01,1];
PARS.max=[1000,50];
PARS.step=0.0001;
PARS.init=[100,5];

%*************Step 4. MHMCMC options *************
%
%The MCMC options can be prescribed in the MCO structure (all optional). These options include
%
%MCO.nout: (>0; default =1000);
%number of parameter vector samples
%
%MCO.silent: (0,1; default = 0);
%%silent code: set to 1 to avoid displaying messages 
%
%MCO.nadapt (>0; default = 50);
%Step-size adaptation frequency: step-size is adjusted based on N samples
%to optimize MHMCMC performance.
%
%MCMC.printfrequency: (>0; default = 10);
%Display MHMCMC progress every N steps 
%
%***Additional MCO options***
%MCO.STOPBEST: (0,1; default 0);
%If set to 1, the MHMCMCM stops if probability = 1; This option should only be used when MLF returns
%0 values (i.e. log probability = 0).
%
%MCO.MH: (0 or 1; default = 1);
%metropolis hastings algorithm (recommended = 1). Setting this option to
%zero results in defficient exploration of parameter space, and incorrect
%estimates of parameter uncertainty.
%
%Example: 
MCO.nout=100000;
MCO.printrate=10000;
MCO.samplerate=20;
MCO.nadapt=2;
MCO.minimumstep=0.001;
MCO.mcmcid=119;

%*************STEP 5. run MHMCMC*************
disp('Optimizing parameters p1 and p2 in the following model, based on observations')
disp('Model: MCH4 = p1 * p2.^(T/10)')
%MHMCMC run
MCMCoutput=MHMCMC(DATA,MLF,PARS,MCO);
%*************STEP 6. Plot results*************
%
%Example: 
PLOT_RESULTS(MCMCoutput,DATA)

%Contact:
%Anthony Bloom: abloom@jpl.nasa.gov
%
%References:
%Bloom, A. A., and M. Williams. "Constraining ecosystem carbon dynamics in a data-limited world: integrating ecological" common sense" in a model?data fusion framework." Biogeosciences 12.5 (2015): 1299-1315.
%Ziehn, T., Marko Scholze, and W. Knorr. "On the capability of Monte Carlo and adjoint inversion techniques to derive posterior parameter uncertainties in terrestrial ecosystem models." Global Biogeochemical Cycles 26.3 (2012).
%Xu, Tao, et al. "Probabilistic inversion of a terrestrial ecosystem model: Analysis of uncertainty in parameter estimation and model prediction." Global Biogeochemical Cycles 20.2 (2006).

end

function PLOT_RESULTS(MCMCoutput,DATA)
%This function plots the MHMCMC outputs.
%Note: Only the final 50% of the samples are used for plotting. From these,
%only 1 in 10 of the samples is used (Ziehn et al., 2012) to minimize the
%correlation between samples.
%
%Top-left panel: parameter 1 probability density function
%Top-right panel: parameter 1 probability density function
%Bottom-left: Model and observations
%Bottom-right: Probability and parameter samples


%Step 1. Keep latter 50% of samples, and sub-sample (e.g. 1 in 10).
%
%All parameter samples
PARSOUT=MCMCoutput.PARSOUT(end/2+1:end,:);
%All associated probabilities
PROB=MCMCoutput.PROB(end/2+1:end);



%Parameter PDFs
figure(1);clf
%"trace plot" for p1, p2 and probability of parameter vector  of each sample
subplot(6,2,8);set(gca,'FontSize',7);plot(PARSOUT(:,1));ylim([0,20]);ylabel('p1');
subplot(6,2,10);set(gca,'FontSize',7);plot(PARSOUT(:,2));ylim([0,10]);ylabel('p2');
subplot(6,2,12);set(gca,'FontSize',7);plot(PROB);ylim([-10,0]);ylabel('Prob.');
xlabel('Sample');


%MHMCMC only returns parameter samples
%To obtain model output as a function of these parameters,
%re-run model with parameter vectors.
%
for n=1:numel(PROB)
    MCH4(n,:)=WETLAND_CH4_MODEL(DATA,PARSOUT(n,:));
end
%plot model median and inter-quartile range
subplot(2,2,3);set(gca,'FontSize',7);
hold on
%model median (based on all model outputs)
h(1)=plot(1:12,median(MCH4,1),'b','LineWidth',2);
%inter-quartile range
h(2)=plot(1:12,prctile(MCH4,25),'b--');plot(1:12,prctile(MCH4,75),'b--');
%Plot observations
h(3)=plot(1:12,DATA.obs,'sr','MarkerSize',5);
%Plot observation uncertainty
for n=1:12;h(4)=plot([n,n],DATA.obs(n)+DATA.obsunc(n).*[1,-1],'r-');end
%legend, labels, month axis, etc.
legend(h,'Model (median)','Model (25th,75th %ile)','Obs','Obs unc.','Location','NorthWest');
ylabel('CH4 flux [mg/m2/day]');xlabel('Month')
set(gca,'xtick',1:12,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
ylim([10,55])

%Making histograms based on parameter (p1 and p2) samples
%Parameters probability density function (based on latter-half of parameter samples)
%
%Note: 
%- log-transform (log10()) and corresponding re-labeling of xticks in code below is only for plotting & visualization purposes
%- Actual parameter values are PARSOUT(:,1) and PARSOUT(:,2) 
%- Log-transform step can be very useful for vizualizing parameters spannign one or more orders of magnitude
%
%
%Parameter 1
subplot(2,2,1);set(gca,'FontSize',7);
hist(log10(PARSOUT(:,1)))
set(gca,'xticklabel',round(10.^(get(gca,'xtick'))*10)/10)
xlabel('Parameter 1 (p1)'); ylabel('Probability density')
%Parameter 2
set(gca,'ytick',[])
subplot(2,2,2);set(gca,'FontSize',7);hist(log10(PARSOUT(:,2)-1))
set(gca,'xticklabel',round((10.^(get(gca,'xtick'))+1)*10)/10)
xlabel('Parameter 2 (p2)');ylabel('Probability density')
set(gca,'ytick',[]);





end



%Model-likelihood function
function P=MODEL_LIKELIHOOD_FUNCTION(DATA,pars)
%This function runs the model using "pars" and compares model
%output with observations. The function returns the probability of "pars".
%Step 1. model run
MCH4 = WETLAND_CH4_MODEL(DATA,pars);

%Step 2. compare to observations
%Note: P is the log probability
%The following is based on the likelihood function
%P = exp(-0.5*(sum((Model - Obs).^2./Unc.^2)));
%See Bloom et al., 2015 and Ziehn et al., 2012 (and references therein) for
%details.
%Correction made on Sep 3, 2020
P = -0.5*(sum((MCH4-DATA.obs).^2./DATA.obsunc.^2));
end

%Model function
function MCH4=WETLAND_CH4_MODEL(DATA,pars)
%Model run, as a function of parameter vector "pars"
MCH4=pars(1)*pars(2).^(DATA.temp/10);
end

