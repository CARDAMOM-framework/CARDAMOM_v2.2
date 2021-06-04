function CBF=cardamomfun_clear_cbf_obs(CBF)

%Overriding obs
of=fieldnames(CBF.OBS);
for n=1:numel(of);CBF.OBS.(of{n})=[];end

%Overrinfding obs unc
of=fieldnames(CBF.OBSUNC);
for n=1:numel(of);
    of2=fieldnames(CBF.OBSUNC.(of{n}));
    for nn=1:numel(of2);
        if isnumeric(CBF.OBSUNC.(of{n}).(of2{nn}));
    CBF.OBSUNC.(of{n}).(of2{nn})=-9999;end
    end
end



%Overrinfding other obs
of=fieldnames(CBF.OTHER_OBS);
for n=1:numel(of);
    of2=fieldnames(CBF.OTHER_OBS.(of{n}));
    for nn=1:numel(of2);
        if isnumeric(CBF.OTHER_OBS.(of{n}).(of2{nn}));
    CBF.OTHER_OBS.(of{n}).(of2{nn})=-9999;end
    end
end


CBF.PARPRIORS=CBF.PARPRIORS*0-9999;
CBF.PARPRIORUNC=CBF.PARPRIORUNC*0-9999;


end