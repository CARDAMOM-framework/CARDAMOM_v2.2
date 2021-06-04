function CBROUT=cardamomfun_combine_parameter_chains(CBR)
%CBR is the standard output from CARDAMOM_RUN_MDF
%
%Last shared on Nov 27, 2019

if isfield(CBR,'chainid');disp('Note: deleting existing "chainid" field');    CBR=rmfield(CBR,'chainid');end

fieldnames=fields(CBR(1,1));
CBROUT=CBR(1);
CBROUT.chainid=CBR(1).PROB*0+1;
    for n=2:numel(CBR)
        for f=1:numel(fieldnames)
        CBROUT.(fieldnames{f})=[CBROUT.(fieldnames{f});CBR(n).(fieldnames{f})];            
        end
    CBROUT.chainid=[CBROUT.chainid;CBR(n).PROB*0+n];
end
    

% 
% PARS=[];
% for n=1:numel(CBR)
% PARS=[PARS;CBR(n).PARS(end/2+1:end,:)];
% end




end