function RT=cardamomfun_pars2RT(CBF,CBR);
disp('Under construction, currently works for C pool purposes...');
%Load inputs and outputs 
%MODEL ID = 4
[FLUXinout]=DALECFUN_FLUX_POOL_IN_OUT(CBF.ID);
%Live biomass inputs
LBin=sum(FLUXinout(:,1:4),2)==1;
%DOM inputs
DOMin=sum(FLUXinout(:,5:6),2)==1;
%Total biomass inputs
TBin=sum(FLUXinout(:,1:6),2)==1;
 

Nyrs=size(CBR.GPP,2)/12;
 
%Single pool RT
for p=1:size(CBR.POOLS,3)
    %Daily means 
    OUTPUTS(:,p)=mean(sum(CBR.FLUXES(:,:,FLUXinout(:,p)==-1),3),2);
    %PGinout(:,p)=(CBR.POOLS(:,end,p)-CBR.PARS(:,17+p))/(365.25*Nyrs);
    %PRT(:,p)=mean(CBR.POOLS(:,:,p),2)./(INPUTS(:,p)-PGinout(:,p))/365.25;    
    PRT(:,p)=mean(CBR.POOLS(:,:,p),2)./(OUTPUTS(:,p))/365.25;    
end
 



 
 
%Store aggregate pool RTs:  live biomass, DOM and total residence times
%Live biomass
AP_INPUTS(:,1)=mean(sum(CBR.FLUXES(:,:,LBin),3),2);
AP_PGinout(:,1)=(sum(CBR.POOLS(:,end,1:4),3)-sum(CBR.PARS(:,18:21),2))/(365.25*Nyrs);
AP_PRT(:,1)=mean(sum(CBR.POOLS(:,:,1:4),3),2)./(AP_INPUTS(:,1)-AP_PGinout(:,1))/365.25;    
%Compare to standard method (done!)
AP_INPUTS(:,2)=mean(sum(CBR.FLUXES(:,:,DOMin),3),2);
AP_PGinout(:,2)=(sum(CBR.POOLS(:,end,5:6),3)-sum(CBR.PARS(:,22:23),2))/(365.25*Nyrs);
AP_PRT(:,2)=mean(sum(CBR.POOLS(:,:,5:6),3),2)./(AP_INPUTS(:,2)-AP_PGinout(:,2))/365.25;    
%Total biomass
AP_INPUTS(:,3)=mean(sum(CBR.FLUXES(:,:,TBin),3),2);
AP_PGinout(:,3)=(sum(CBR.POOLS(:,end,1:6),3)-sum(CBR.PARS(:,18:23),2))/(365.25*Nyrs);
AP_PRT(:,3)=mean(sum(CBR.POOLS(:,:,1:6),3),2)./(AP_INPUTS(:,3)-AP_PGinout(:,3))/365.25;

RT.LB=AP_PRT(:,1);
RT.DOM=AP_PRT(:,2);
RT.ECO=AP_PRT(:,3);

RT.C_POOLS=PRT;


end