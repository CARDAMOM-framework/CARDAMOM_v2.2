function cardamom_vvuq_cbf_summary(CBF)
%function cardamom_vvuq_cbf_summary(CBF)
%Summarizes contents of CBF file
%
%Written by A Bloom, Thu Apr 11, 2019

%Switch is in place in case later version met fields are not compatible!
figure(1);clf




MD=CARDAMOM_MODEL_LIBRARY(CBF.ID);

for s=1:9;subplot(3,3,s);set(gca,'FontSize',6);end

subplot(3,3,1);plot(CBF.MET(:,1));
ylabel('Time [Days since Jan 01, 2001]');xlabel('Timestep (months)')
subplot(3,3,2);plot(CBF.MET(:,2));
ylabel('min temp [C]');xlabel('Timestep (months)')
subplot(3,3,3);plot(CBF.MET(:,3));
ylabel('max temp [C]');xlabel('Timestep (months)')
subplot(3,3,4);plot(CBF.MET(:,4));
ylabel('Solar irradiance [MJ/m2/day]');xlabel('Timestep (months)')
subplot(3,3,5);plot(CBF.MET(:,5));
ylabel('CO2');xlabel('Timestep (months)')
subplot(3,3,6);plot(CBF.MET(:,6));
if MD.nomet>6
ylabel('Day of year');xlabel('Timestep (months)')
subplot(3,3,7);plot(CBF.MET(:,7));
if MD.nomet>7
ylabel('Burned area [m2/m2]');xlabel('Timestep (months)')
subplot(3,3,8);plot(CBF.MET(:,8));
ylabel('VPD [hPa]');xlabel('Timestep (months)')
subplot(3,3,9);plot(CBF.MET(:,9));
ylabel('Precip. [mm/day]');xlabel('Timestep (months)');
end
end






figure(2);clf
ofn=fields(CBF.OBS);
for n=1:numel(ofn)
    subplot(3,4,n);
    
    %plotting 
    
    plot(fill2nan(CBF.OBS.(ofn{n})));title(ofn{n});
    if isempty(CBF.OBS.(ofn{n}));set(gca,'Color','None');end
end

disp('*********OBS UNC*********')
oufn=fields(CBF.OBSUNC);
for n=1:numel(oufn)
disp(sprintf('***%s***',oufn{n}))
fn=fields(CBF.OBSUNC.(oufn{n}));
for nn=1:numel(fn)
    if isstr(CBF.OBSUNC.(oufn{n}).(fn{nn}))==0;
disp(sprintf('%s =%g',strjust(sprintf('%50s',fn{nn})),CBF.OBSUNC.(oufn{n}).(fn{nn})));
    end
end
end
    

disp('*********OTHER OBS*********')
oufn=fields(CBF.OTHER_OBS);
for n=1:numel(oufn)
disp(sprintf('***%s***',oufn{n}))
fn=fields(CBF.OTHER_OBS.(oufn{n}));
for nn=1:numel(fn)
    if isstr(CBF.OTHER_OBS.(oufn{n}).(fn{nn}))==0;
        disp(sprintf('%s =%g',strjust(sprintf('%30s',fn{nn})),CBF.OTHER_OBS.(oufn{n}).(fn{nn})));
    end
end
end

%list parameters
MD=CARDAMOM_MODEL_LIBRARY(CBF.ID);
disp('Parameter priors')
for n=1:MD.nopars
    if CBF.PARPRIORS(n)>-9999 | CBF.PARPRIORUNC(n)>-9999
        disp(sprintf('%s: Prior = %2.2f; Unc = %2.2f',MD.parname{n},CBF.PARPRIORS(n),CBF.PARPRIORUNC(n)));
    end
end

end


function [v,pts]=fill2nan(v)

pts=find(v~=-9999);
v(v==-9999)=NaN;

end