function [FLUXinout]=DALECFUN_FLUX_POOL_IN_OUT(ID)
%This function loads two binary arrays with flux inputs/outputs for each
%pool

if ID==4;
    FLUXinout=dalec2fires_FLUXinout;
elseif ID==805;
    FLUXinout=dalec805_FLUXinout;
    elseif ID==809;
    FLUXinout=dalec805_FLUXinout;
        elseif ID==811;
    FLUXinout=dalec805_FLUXinout;
    
end
    



end

function FLUXinout=dalec2fires_FLUXinout

% DR.GPP=permute(FLUXES(1,:,:),[3,2,1]); 
% DR.temprate=FLUXES(2,:);
% DR.respiration_auto=FLUXES(3,:); PIN
% DR.leaf_production=FLUXES(4,:); PIN(4,:)=[0,1,0,0,0,0];
% DR.labile_production=FLUXES(5,:); PIN(
% DR.root_production=FLUXES(6,:);
% DR.wood_production=FLUXES(7,:);
% DR.labile_release=FLUXES(8,:);
% DR.leaffall_factor=FLUXES(9,:);
% DR.leaflitter_production=FLUXES(10,:);
% DR.woodlitter_production=FLUXES(11,:);
% DR.rootlitter_production=FLUXES(12,:);
% DR.respiration_het_litter=FLUXES(13,:);
% DR.respiration_het_som=FLUXES(14,:);
% DR.litter2som=FLUXES(15,:);
% DR.labrelease_factor=FLUXES(16,:);
% DR.fires=FLUXES(17,:);
% DR.fire_per_pool = FLUXES(18:23); %pool to atmosphere
% DR.fire_transfers = FLUXES(24:28); % pool to litter (and litter to SOM).
FLUXinout(1,:)=[0,0,0,0,0,0];%GPP (total)
FLUXinout(2,:)=[0,0,0,0,0,0];%temp
FLUXinout(3,:)=[0,0,0,0,0,0];%auto
FLUXinout(4,:)=[0,1,0,0,0,0];%fol prod
FLUXinout(5,:)=[1,0,0,0,0,0];%lab prod
FLUXinout(6,:)=[0,0,1,0,0,0];%root prod
FLUXinout(7,:)=[0,0,0,1,0,0];%wood prod
FLUXinout(8,:)=[-1,1,0,0,0,0];%lab release
FLUXinout(9,:)=[0,0,0,0,0,0];%leaffall factor
FLUXinout(10,:)=[0,-1,0,0,1,0];%leaflitter
FLUXinout(11,:)=[0,0,0,-1,0,1];%woodlitter
FLUXinout(12,:)=[0,0,-1,0,1,0];%rootlitter
FLUXinout(13,:)=[0,0,0,0,-1,0];%litter resp
FLUXinout(14,:)=[0,0,0,0,0,-1];%soil resp
FLUXinout(15,:)=[0,0,0,0,-1,1];%litt2som
FLUXinout(16,:)=[0,0,0,0,0,0];%labrelease factor
FLUXinout(17,:)=[0,0,0,0,0,0];%fires (total)
FLUXinout(18,:)=[-1,0,0,0,0,0];%labburned
FLUXinout(19,:)=[0,-1,0,0,0,0];%GPP
FLUXinout(20,:)=[0,0,-1,0,0,0];%GPP
FLUXinout(21,:)=[0,0,0,-1,0,0];%GPP
FLUXinout(22,:)=[0,0,0,0,-1,0];%GPP
FLUXinout(23,:)=[0,0,0,0,0,-1];%GPP
FLUXinout(24,:)=[-1,0,0,0,1,0];%GPP
FLUXinout(25,:)=[0,-1,0,0,1,0];%GPP
FLUXinout(26,:)=[0,0,-1,0,1,0];%GPP
FLUXinout(27,:)=[0,0,0,-1,0,1];%GPP
FLUXinout(28,:)=[0,0,0,0,-1,1];%GPP






end


function FLUXinout=dalec805_FLUXinout

% DR.GPP=permute(FLUXES(1,:,:),[3,2,1]); 
% DR.temprate=FLUXES(2,:);
% DR.respiration_auto=FLUXES(3,:); PIN
% DR.leaf_production=FLUXES(4,:); PIN(4,:)=[0,1,0,0,0,0];
% DR.labile_production=FLUXES(5,:); PIN(
% DR.root_production=FLUXES(6,:);
% DR.wood_production=FLUXES(7,:);
% DR.labile_release=FLUXES(8,:);
% DR.leaffall_factor=FLUXES(9,:);
% DR.leaflitter_production=FLUXES(10,:);
% DR.woodlitter_production=FLUXES(11,:);
% DR.rootlitter_production=FLUXES(12,:);
% DR.respiration_het_litter=FLUXES(13,:);
% DR.respiration_het_som=FLUXES(14,:);
% DR.litter2som=FLUXES(15,:);
% DR.labrelease_factor=FLUXES(16,:);
% DR.fires=FLUXES(17,:);
% DR.fire_per_pool = FLUXES(18:23); %pool to atmosphere
% DR.fire_transfers = FLUXES(24:28); % pool to litter (and litter to SOM).
FLUXinout(1,:)=[0,0,0,0,0,0,0];%GPP (total)
FLUXinout(2,:)=[0,0,0,0,0,0,0];%temp
FLUXinout(3,:)=[0,0,0,0,0,0,0];%auto
FLUXinout(4,:)=[0,1,0,0,0,0,0];%fol prod
FLUXinout(5,:)=[1,0,0,0,0,0,0];%lab prod
FLUXinout(6,:)=[0,0,1,0,0,0,0];%root prod
FLUXinout(7,:)=[0,0,0,1,0,0,0];%wood prod
FLUXinout(8,:)=[-1,1,0,0,0,0,0];%lab release
FLUXinout(9,:)=[0,0,0,0,0,0,0];%leaffall factor
FLUXinout(10,:)=[0,-1,0,0,1,0,0];%leaflitter
FLUXinout(11,:)=[0,0,0,-1,0,1,0];%woodlitter
FLUXinout(12,:)=[0,0,-1,0,1,0,0];%rootlitter
FLUXinout(13,:)=[0,0,0,0,-1,0,0];%litter resp
FLUXinout(14,:)=[0,0,0,0,0,-1,0];%soil resp
FLUXinout(15,:)=[0,0,0,0,-1,1,0];%litt2som
FLUXinout(16,:)=[0,0,0,0,0,0,0];%labrelease factor
FLUXinout(17,:)=[0,0,0,0,0,0,0];%fires (total)
FLUXinout(18,:)=[-1,0,0,0,0,0,0];%labburned
FLUXinout(19,:)=[0,-1,0,0,0,0,0];%GPP
FLUXinout(20,:)=[0,0,-1,0,0,0,0];%GPP
FLUXinout(21,:)=[0,0,0,-1,0,0,0];%GPP
FLUXinout(22,:)=[0,0,0,0,-1,0,0];%GPP
FLUXinout(23,:)=[0,0,0,0,0,-1,0];%GPP
FLUXinout(24,:)=[-1,0,0,0,1,0,0];%GPP
FLUXinout(25,:)=[0,-1,0,0,1,0,0];%GPP
FLUXinout(26,:)=[0,0,-1,0,1,0,0];%GPP
FLUXinout(27,:)=[0,0,0,-1,0,1,0];%GPP
FLUXinout(28,:)=[0,0,0,0,-1,1,0];%GPP
FLUXinout(29,:)=[0,0,0,0,0,0,-1];%GPP
FLUXinout(30,:)=[0,0,0,0,0,0,-1];%GPP






end