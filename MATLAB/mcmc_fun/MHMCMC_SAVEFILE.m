function MCMC=MHMCMC_SAVEFILE(DATA,MLF,PARS,MCO)
disp('UNTESTED...')
%Based on haario et al. 2001 algorithm
%
%The following adaptive Metropolis-Hasting Markov CHain Monte Carlo (MHMCMC) function
%can be used to sample the probability distribution of a
%user-defined model parameter set. A working example is provided in MHMCMC_EXAMPLE.m
%
%INPUTS:
%
%***DATA: Data structure***
%DATA: all fields are user defined; all model datasets, observations, and
%uncertainties should be stored in the DATA structure
%
%***MLF: model likelihood function***
%must accept DATA and pars vector, and must return the log likelihood, i.e.
%P=MLF(DATA,pars); P is the log likelihood.
%
%***PARS: Parameter info structure***
%Contains minimum and maximum parameter values: PARS.min,
%PARS.max. PARS.step has a default value of 0.001 for each parameter;
%Fields:
%PARS.init = initial values
%PARS.min = minimum values
%PARS.max = maximum values
%PARS.step = initial step size
%
%***MCO: MCMC  optionsoptions.***
%MCO.nout: number of parameter vector samples (>0; default =1000);
%MCO.silent: set to 1 to avoid displaying messages  (0,1; default = 0);
%MCO.nadapt: Step-size adaptation frequency. (>0; default = 50);
%MCMC.printrate: Display MHMCMC progress every N steps (1,inf; default = 10);
%MCO.STOPBEST: If set to 1, the MHMCMCM stops if probability = 1; (0,1; default 0);
%MCO.MH: Switch MHMCMC on (1) or off (0); default = 1.
%
%Outputs:
%
%***MCMC: Output structure***
%Contains output structure including parameter samples (PARSOUT),
%probability of parameter samples (PROB).
%
%Function last edited on Mar 18th 2016 by A. Anthony Bloom
%abloom@jpl.nasa.gov
%
%NOTE: ongoing edits: prescribing PCA rotation on step distribution





%Default MCOpt 
if nargin<4; MCO=[];end
%silent code
if isfield(MCO,'savefile')==0;MCO.savefile=[];end

if isempty(MCO.savefile) | isfile(MCO.savefile)==0;

if isfield(MCO,'silent')==0;MCO.silent=0;end
%metropolis hastings algorithm
if isfield(MCO,'MH')==0;MCO.MH=1;end
%number of samples
if isfield(MCO,'nout')==0;MCO.nout=1000;end
%Adapt every so often (memory limited)
if isfield(MCO,'nadapt')==0;MCO.nadapt=ceil(MCO.nout/20000);end

%Save every 100 adaptations
if isfield(MCO,'saverate')==0;MCO.saverate=MCO.nadapt*100;end
%sample frequency
if isfield(MCO,'samplerate')==0;MCO.samplerate=ceil(MCO.nout/2000);end
%minimum step size: only adjust for large number of parameters
if isfield(MCO,'minstepsize')==0;MCO.minstepsize=0.001;end


%Stop if probability of 1 is found
if isfield(MCO,'STOPBEST')==0;MCO.STOPBEST=0;end
%Display every N adaptations
if isfield(MCO,'printrate')==0;MCO.printrate=10;end
%Graphical display of MCMC
if isfield(MCO,'graphical')==0;MCO.graphical=0;end



%If no initial parameter vector (PARS.init) was prescribed, random
%parameters within prescribed parameter range are chosen
if isfield(PARS,'init')==0;disp('No initial values provided, randomly selecting starting point');PARS.init=norm2pars(PARS,rand(size(PARS.min)));end

    
%PCA rotation step 
    

%displaying initial parameter values
if MCO.silent==0;    
    disp('Initial parameter values found');
   % disp(sprintf('%3d %10.5f\n ',[1:length(PARS.init);PARS.init]));
end







%pars0 and P0
pars0=PARS.init;
%Test Model-Likelihood function with inital parameter vector pars0
P0=MLF(DATA,pars0);
%Warning in case initial parameter vector pars0 yields probability of 0.
if isfinite(P0)==0 | P0>0;disp(sprintf('Warning - initial probability (P0) = %f',P0));end


Nsamples=MCO.nout/MCO.nadapt;
Nsamplesout=floor(MCO.nout/MCO.samplerate);
Npars=numel(pars0);

%Step 
STEP.covariance=zeros(Npars);
STEP.cholesky=zeros(Npars);

%declaring local parameter sample array
PARSALL=repmat(pars0*0,[Nsamples,1]);
%declaring full sample (subsampled) parameter array
PARSout=repmat(pars0*0,[Nsamplesout,1]);
%probability
PROBout=zeros(1,Nsamplesout)*NaN;

%loop counters:
%declaring (number of samples), local acceptance, and number of raw runs
n=0;
local_acc=0;
runcount=0;


STINDEX=1;
else 
    
load(MCO.savefile);

    STINDEX=n+1;
    
end


for n=STINDEX:MCO.nout;

    %new parameter sample
    %pstep=[PARS.PCArotation*[PARS.step.*randn(size(pars0))]']';
    rv1=randn(1,Npars);rv2=randn(1,Npars);
    pstep=rv1*STEP.cholesky + MCO.minstepsize*rv2;
    %pars=norm2pars(PARS,pars2norm(PARS,pars0)+pstep);
    pars=norm2pars(PARS,pars2norm(PARS,pars0)+pstep);
    
    %probability is now log likelihood
    %NOTE: probability is zero if parameters are not within prior range
    P=MLF(DATA,pars)+log(double(sum(PARS.min>pars)==0))+log(double(sum(PARS.max<pars)==0));
    %P and P0 checks
% 
%     if MCO.graphical==1
%         %plot 2 dims only
%         plot3(log10([pars0(1),pars(1)]),log10([pars0(2),pars(2)]),[0,0],'.--','LineWidth',0.5,'Color',ones(1,3)*0.8);
%         drawnow;
%     end
    
   %Accept-Reject based on probability ratio and random number
   %(see Ziehn et al., 2012)
    if P-P0>log(rand(1))^MCO.MH
         %local loop constant
        local_acc=local_acc+1;
        
        %pars0
        pars0=pars;
        %As of March 18, 2016: probability is no longer a standard output
        %of the MHMCMC function
        P0=P; 
        %Storing in PARSOUT 

    end
    
    %Storing outputs
        if mod(n,MCO.samplerate)==0
        Nssamples=n/MCO.samplerate;

        PARSout(Nssamples,:)=pars0;
    PROBout(Nssamples)=P0;

        end
        
    if mod(n,MCO.nadapt)==0
        Nssamples=n/MCO.nadapt;
    %Accumulate parameters
    PARSALL(Nssamples,:)=pars0;
  

    %Adjust step size, re-setting acceptance count and displaying
    %performance

    if  MCO.silent==0  & mod(n,MCO.printrate)==0;
        disp(['Acceptance rate = ',num2str(local_acc/MCO.nadapt*100),'%']);
        disp(['Accepted ',num2str(n),' solutions'])
        disp(['(Log) Probability = ',num2str(P0)])
        disp(['mean step size = ',num2str(sqrt(mean(diag(STEP.covariance))))])
        if Nssamples<10*Npars;   disp(sprintf('%2.2f%% of samples aquired for covariance sampler...',Nssamples/10/Npars*100));end
    end
        
    
        %Adjusting step size for optimal MHMCMC performance
        inccov=1;
        if n<MCO.nout/2  & Nssamples>10*Npars;     
            %If cov is empty
            if STEP.covariance(1,1)==0 | inccov==0
                STEP.covariance=cov(pars2norm(PARS,PARSALL(floor(Nssamples/2)+1:Nssamples,:)));
                STEP.runmean=mean(pars2norm(PARS,PARSALL(floor(Nssamples/2)+1:Nssamples,:)));
                STEP.covsamples=Nssamples-floor(Nssamples/2);
            else
                %Calculate covariance of parameters
                %[CMOUT,Mi]=statfun_add_sample_to_covmat(CM,N,M, x,ar)
                [STEP.covariance,STEP.runmean]=statfun_add_sample_to_covmat_local(STEP.covariance,STEP.covsamples,STEP.runmean,pars2norm(PARS,PARSALL(Nssamples,:)),1);
                STEP.covsamples= STEP.covsamples+1;
                if mod(Nssamples,2)==0;
                [STEP.covariance,STEP.runmean]=statfun_add_sample_to_covmat_local(STEP.covariance,STEP.covsamples,STEP.runmean,pars2norm(PARS,PARSALL(Nssamples/2,:)),-1);
                STEP.covsamples= STEP.covsamples-1;
                end
                
            end
            
            STEP.cholesky=chol(STEP.covariance*2.38^2/Npars);
        end
        %Re-setting acceptance rate count
        local_acc=0;
    end
    
    if isempty(MCO.savefile)==0 & mod(n,MCO.saverate)==0;
    save(MCO.savefile);
    end
    

    %stop loop if option to stop when log probability (P)=0 has been selected
    if (MCO.STOPBEST==1 & P==0)
        disp('Found a P=1 solution!')
        MCO.nout=n;
    end
%End of MHMCMC main loop
end


%MCMC success status




    disp('MHMCMC successfully completed!')
    mcmcsuccess=1;




% [PARSOUT,PROB,PARS,mcmcsuccess]=MHMCMC(DATA,MLF,PARS,MCO)
%Outputs
%Pameter vectors
MCMC.PARSOUT=PARSout;
%Probability of parameter vectors
MCMC.PROB=PROBout;
%Parameter info (same as input + additional fields)
MCMC.PARS=PARS;
%Success status
MCMC.success=mcmcsuccess;
%Options
MCMC.MCO=MCO;




end

%normalize pars w.r.t. minimum and maximum parameter values
function npars=pars2norm(PARS,pars)
for n=1:size(pars,1);
    npars(n,:)=(log(pars(n,:)./PARS.min))./(log(PARS.max./PARS.min));
end
end

%scale normalized parameters to true values w.r.t. minimum and maximum
%parameter values
function pars=norm2pars(PARS,npars)
pars=PARS.min.*(PARS.max./PARS.min).^npars;
end


%Adjust the step size in all dimensions every N iterations to optimize MHMCMC performance
function step=adjuststepsize(PARS,acc,ai,pars,minstepsize)
%The step size "step" is adjusted here according to acceptance stats "acc".
%See Bloom et al., 2015 for details.
%
step=PARS.step;

%Adaptation factors
%adapfac_all=1.5;
%adapfac_step=1.5;
adapfac_all=2;
adapfac_step=1.1;
fac=2;



    if acc>0 & acc<round(0.23*ai) 
         step=step*(1/adapfac_all);
    elseif acc>(0.44*ai)
      step=step*adapfac_all;
    end

    %Step-specific adaptations
    %For PCA comparison, need to either (a) rotate steps to par-space
    %OR, (b) rotate pars to step-space (correct, since this then allows
    %step to be adapted w.r.t. multivariate drift)
    
    %PARS.PCArotation: reverse rotation
    if acc>3;stdacc=std( pars2norm(PARS,pars(1:end,:))*PARS.PCArotation);end
    
for n=1:length(step)


    
    %Adapting step to mitigate "drift"
    if acc>3
    if step(n)>stdacc(n);
    step(n)=step(n);%/1.1;
    elseif step(n)<stdacc(n)/fac & acc<(0.23*ai);
      step(n)=step(n)*adapfac_step;%minstepsize;
    end
    end
    
        %ensuring step is not too big or small
    if step(n)>1
         step(n)=step(n)*(1/adapfac_step);
    elseif step(n)<minstepsize
      step(n)=minstepsize;%step(n)*adapfac_step;
    end
    
end

step(step>1)=1;

end

function [CMOUT,Mi]=statfun_add_sample_to_covmat_local(CM,N,M, x,ar)
Mi=(M*N+ ar*x)/(N+ar);
CMOUT=CM*(N-1)/(N-1+ar)+((N)* M'*M- (N+ar)*Mi'*Mi + x'*x*ar)/(N-1+ar);
end
