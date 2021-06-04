function gr=cardamomfun_convergence_tests(A1,A2,A3)
%

%CASE 1: code recognizes CBF and CBR structure
%In which case inputs are (CBF,CBR,chains), where "chains" is a vector
%denoting which chains to intercompare (e.g. intercomparison between all
%four chains is [1,1,1,1], intecomparison between first and third only is
%[1,0,1,0]
%
%Example use:
% MCO.niterations=1e7;
% for n=1:3;CBR(n)=CARDAMOM_RUN_MDF(CBF,MCO);end
% CBR=cardamomfun_combine_parameter_chains(CBR);
% [gr]=cardamomfun_convergence_tests(CBF,CBR,[1,1,1])
%
%Last updated by Anthony Bloom (07/10/2020)

if isstruct(A1) & isfield(A1,'MET')==1
    CBF=A1;
    CBR=A2;
    chains=A3;
    MD=CARDAMOM_MODEL_LIBRARY(CBF.ID);

[gr]=cardamomfun_convergence_tests_cbf_cbr_structs(MD,CBR.PARS,CBR.chainid,chains,1);

%CASE 2: code recognizes PXI structure, only one argument required
%
elseif isstruct(A1) & isfield(A1,'run_name');
    %Recognizes PXI format[gr]=cardamomfun_convergence_tests
    PXI=A1;
    [gr]=cardamomfun_convergence_tests_pxi(PXI);
    
end
    
    



end

function     PXI=cardamomfun_convergence_tests_pxi(PXI)

%Step 1. loop through all cbr files
    MD=CARDAMOM_MODEL_LIBRARY(PXI.ID);
    
PXI.chain_convergence_max_gr_map=PXI.area*NaN;

for n=1:size(PXI.cbrfilename)
    cbrfiles=PXI.cbrfilename(n,PXI.chains(n,:)>0);
    
    if isempty(cbrfiles)==0
    %open cbr files
    [PARS,ANC]=CARDAMOM_READ_BINARY_FILEFORMAT(cbrfiles,MD);
    chains=zeros(max(unique(ANC.chainid)));
    chains(unique(ANC.chainid))=1;
    if numel(chains)>1
    [gr]=cardamomfun_convergence_tests_cbf_cbr_structs(MD,PARS,ANC.chainid,chains,0);
    else
        gr=nan(1,MD.nopars);
    end

PXI.chain_convergence(n,:)=gr;
        PXI.chain_convergence_max_gr_map(PXI.r(n),PXI.c(n))=max(PXI.chain_convergence(n,:));

    else
        PXI.chain_convergence(n,:)=nan(1,MD.nopars);
    end
    
    
end

%Add map here


% PXI=cardamomfun_convergence_tests_pxi(PXI);



end

function [gr]=cardamomfun_convergence_tests_cbf_cbr_structs(MD,PARS,chainid,chains,plotfig)

gr=[];
if MD.ID~=101;
PARS(:,[12,15])=mod(PARS(:,[12,15]),365.25)+365.25;
else
    PARS(:,[6,9])=mod(PARS(:,[6,9]),365.25)+365.25;
end
    
NPARSall=MD.par2nor(PARS);

k=1;N=size(NPARSall,1);%nc=total(chains);

%total number of elements from all chains
for n=1:max(chainid)
    chainpts{n}=find(chainid==n);
    nelements(n)=numel(chainpts{n});
end
minelements=min(nelements(chains==1));
%

for n=1:numel(chains);
    if chains(n)==1;
        pts=chainpts{n}(1:minelements);
        NPARS(:,:,k)=NPARSall(pts,:);
        k=k+1;
    end
end

for n=1:MD.nopars
gr(n)=mcmc_psrf(NPARS(:,n,:));
end






%Making figures
if plotfig==1;
%Convergence
figure(1);clf
loadcolornames
cols={'r','b','k',forestgreen};
cols=repmat(cols,[1,ceil(numel(chains)/4)]);
for n=1:MD.nopars;subplot(7,6,n);set(gca,'FontSize',4);

    for nn=1:size(NPARS,3);ph=plotpdf(NPARS(:,n,nn),cols{nn}); hold on;        set(ph,'LineWidth',0.5);end
    if max(chainid)>1
text(0.98,1.14,sprintf('(%2.2f)',gr(n)),'Units','Normalized','HorizontalAlignment','right','FontSize',6,'FontWeight','Bold','Color','r');
title(sprintf('P%i   ',n),'FOntSize',6);
set(gca,'yticklabel','');
    end
end


%Convergence
figure(2);clf
cols={'r','b','k',forestgreen};
cols=repmat(cols,[1,ceil(numel(chains)/4)]);
for n=1:MD.nopars;subplot(7,6,n);set(gca,'FontSize',4);hold on

    for nn=1:size(NPARS,3);ph=plot(NPARS(:,n,nn),'Color',cols{nn}); hold on;        set(ph,'LineWidth',0.5);end
    ylim([0,1]);
    if max(chainid)>1
gr(n)=mcmc_psrf(NPARS(:,n,:));
text(0.98,1.14,sprintf('(%2.2f)',gr(n)),'Units','Normalized','HorizontalAlignment','right','FontSize',6,'FontWeight','Bold','Color','r');
title(sprintf('P%i   ',n),'FOntSize',6);
    end
end

end


end