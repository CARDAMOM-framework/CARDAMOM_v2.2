function varargout=CARDAMOM_WRITE_BINARY_FILEFORMAT(CBF,filename)
%function CARDAMOM_WRITE_BINARY_FILEFORMAT(CBF,filename)
%
%INPUTS: 
% - CBF: CARAMOM binary structure
% - "filename": name of file 
%Writes CBF structure to CARDAMOM binary file (.cbf) format
%See CARDAMOM_READ_BINARY_FILEFORMAT.m for obtaining a "CBF" structure
%template
%
%Last modified by A.A. Bloom 2019/07/24

if strcmp(filename(end-2:end),'cbf');[BINARY_DATA,CBF]=write_cbf_file(CBF,filename);if nargout==2;
    varargout={BINARY_DATA,CBF};
elseif nargout==1;
    varargout={BINARY_DATA};
    end
end

if strcmp(filename(end-2:end),'cbr') | strcmp(filename(end-7:end),'cbrSTART');write_cbr_file(CBF,filename);end
    

end



function varargout=write_cbr_file(CBF,filename,INFO)

%Read CBR file
PARSALL=CBF';




    %number of parameter sets 
    fd=fopen(filename,'w');if fd==-1;disp('WARNING: file not opened, expect error!!');end;fwrite(fd,PARSALL(:),'real*8');fclose(fd);


varargout={''};

end




function [BINARY_DATA,CBF]=write_cbf_file(CBF,filename)

SD=zeros(100,1)-9999;
PR=CBF.PARPRIORS;
PRU=CBF.PARPRIORUNC;


OPR=CBF.RAW.OTHERPRIORS;
OPRU=CBF.RAW.OTHERPRIORSUNC;

%Number of timesteps
nodays=size(CBF.MET,1);
%Number of met fields
nomet=size(CBF.MET,2);
%
CBF.OBS=compile_cbf_obs_fields(CBF);
%Number of time-resolved observations
noobs=size(CBF.OBS,2);

CBF.nomet=nomet;
CBF.noobs=noobs;
CBF.nodays=nodays;

% if CBF.nomet~=nomet;warning('CBF.nomet updated for consistency with CBF.MET dimensions');CBF.nomet=nomet;end
% if CBF.noobs~=noobs;warning('CBF.noobs updated for consistency with CBF.OBS dimensions');CBF.noobs=noobs;;end
% if CBF.nodays~=nodays;warning('CBF.nodays updated for consistency with CBF.MET and CBF.OBS dimensions');CBF.nodays=nodays;end


%Check if nomet is supported by library
MD=CARDAMOM_MODEL_LIBRARY(CBF.ID);
if MD.nomet~=CBF.nomet;error('Wrong number of CBF.MET columns, check MD=CARDAMOM_MODEL_LIBRARY(CBF.ID) for correct number');end


SD(1:7)=[CBF.ID;CBF.LAT;CBF.nodays;CBF.nomet;CBF.noobs;CBF.EDC;CBF.EDCDIAG];
% %EDCDIAG - set "on" by default
% if isfield(MD,'EDCDIAG')==0;SD(7)=1;end

%Write time-invariant terms
[SD,OPRU]=write_obs_uncertainty_fields(CBF,SD,OPRU);

%Write other obs constraints
[OPR,OPRU]=write_other_obs_constraints(CBF,OPR,OPRU);

%Prescribe reference met:
if isfield(CBF,'mmet') & numel(CBF.mmet)==CBF.nomet
    MTEMPDATA=[CBF.mmet,zeros(1,CBF.noobs)-9999];
else 
    MTEMPDATA=[];
end
    
    
ALLTEMPDATA=[[CBF.MET,CBF.OBS];MTEMPDATA]';
%Order = timestep1 (met, obs), timestep2 (met, obs)....
TEMPDATA=ALLTEMPDATA(1:end)';
if prod(size(ALLTEMPDATA))~=(CBF.nomet+CBF.noobs)*(CBF.nodays+double(numel(MTEMPDATA)>0));
    disp('Error!! incorrectly set OBS or MET data vector!!')
else

BINARY_DATA=[SD;PR;PRU;OPR;OPRU;TEMPDATA];

fd=fopen(filename,'w');
if fd==-1;disp('could not create file! soz!');disp(filename);keyboard;end
fwrite(fd,BINARY_DATA,'real*8');
fclose(fd);


end


end



%Time-resolved observations
function CBFOBS=compile_cbf_obs_fields(CBF)


if isfield(CBF,'OBS');

if isnumeric(CBF.OBS); CBFOBS=CBF.OBS;
warning('CARDAMOM_WRITE_BINARY_FILEFORMAT: CBF.OBS matrix used to define observations (will soon be obsolete)');

else


obsnames={'GPP','LAI','NBE','ABGB','ET','EWT','BAND1','BAND2','BAND3','BAND4','SOM','CH4'};
fnames=fieldnames(CBF.OBS);

    for n=1:numel(obsnames);isobs(n)=sum(strcmp(obsnames{n},fnames));end

    nobs=find(isobs==1,1,'last');
        %Defining CBFOBS
        N=size(CBF.MET,1);
    CBFOBS=zeros(N,max(nobs,3))-9999;
for n=find(isobs==1);
    %Empty fields are tolerated (consequently no writing to file)
    if isempty(CBF.OBS.(obsnames{n}))==0 & numel(CBF.OBS.(obsnames{n}))==N;CBFOBS(:,n)=CBF.OBS.(obsnames{n});end
    if numel(CBF.OBS.(obsnames{n}))~=N & isempty(CBF.OBS.(obsnames{n}))==0;warning('CBF.OBS field dimensions not compatible with CBF.MET, writing -9999');keyboard;end
end



end

else
    %FIlling in OBS with -9999 values;  
     %Eventually step can be made obsolete when C code can accept
     %(1) noobs, and (2) obsid
     
        CBFOBS=zeros(size(CBF.MET,1),3)-9999;
end



end

%Time-resolved observation uncertainties
function [SD,OPRU]=write_obs_uncertainty_fields(CBF,SD,OPRU)

%From CARDAMOM_READ_BINARY_DATA.c 
% DATA->nee_annual_unc=statdat[13];
% DATA->et_annual_unc=statdat[14];
% DATA->nee_obs_unc=statdat[15];if (statdat[15]<0){DATA->nee_obs_unc=0.5;}
% DATA->et_obs_unc=statdat[16];if (statdat[16]<0){DATA->et_obs_unc=2;}
% DATA->ewt_annual_unc=statdat[17];
% DATA->ewt_obs_unc=statdat[18];if (statdat[18]<0){DATA->ewt_obs_unc=50;}
% DATA->gpp_annual_unc=statdat[19];
% DATA->gpp_obs_unc=statdat[20];if (statdat[18]<0){DATA->gpp_obs_unc=2;}
% DATA->ch4_annual_unc=statdat[24];  /*shuang*/
% DATA->ch4_obs_unc=statdat[25];if (statdat[25]<0){DATA->ch4_obs_unc=0.5;}  /*shuang*/
% DATA->ch4_obs_threshold=statdat[26]; if (statdat[26]<0){DATA->ch4_obs_threshold=0;}  /*shuang*/

%NBE
SD(14)=CBF.OBSUNC.NBE.annual_unc;
SD(16)=CBF.OBSUNC.NBE.seasonal_unc;
%CBF.OBSUNC.NBE.info='Single point (default = 0.5) and annual (annual_unc) NBE uncertainty [gC/m2/d]';

%ET 
SD(15)=CBF.OBSUNC.ET.annual_unc;
SD(17)=CBF.OBSUNC.ET.unc;
%CBF.OBSUNC.ET.info='Single point (default = 2) and annual (annual_unc) ET uncertainty [mm/d]';

%EWT
SD(18)=CBF.OBSUNC.EWT.annual_unc;
SD(19)=CBF.OBSUNC.EWT.unc;
%CBF.OBSUNC.EWT.info=sprintf('Single point (default = 50) and annual (N/A yet) EWT uncertainty [mm]');

%GPP 
SD(20)=CBF.OBSUNC.GPP.annual_unc;
SD(21)=CBF.OBSUNC.GPP.unc;
%CBF.OBSUNC.GPP.info='Single point (default = 2) and annual (annual_unc) GPP uncertainty [gC/m2/d]';
%gpp abs: option on treating GPP constraint as absolute or relative
SD(8)=CBF.OBSUNC.GPP.gppabs;
%CBF.OBSUNC.GPP.gppabs_info='Set to "1" for GPP data, set to "0" for SIF data"';

%Threshold values (for positive-definite values)
SD(22)=CBF.OBSUNC.ET.obs_unc_threshold;
SD(23)=CBF.OBSUNC.GPP.obs_unc_threshold;

%CH4
SD(25)=CBF.OBSUNC.CH4.annual_unc;
SD(26)=CBF.OBSUNC.CH4.unc;
SD(27)=CBF.OBSUNC.CH4.obs_unc_threshold;

%Time-resolved biomass uncertaiknty
OPRU(2)=CBF.OBSUNC.ABGB.unc;
%CBF.OBSUNC.ABGB.info=sprintf('Uncertainty on time-resolved biomass CBF.OBS.ABGB');
%Time-resolved SOM uncertainty
OPRU(8)=CBF.OBSUNC.SOM.unc;
%CBF.OBSUNC.SOM.info=sprintf('Single point uncertainty on time-resolved SOM, CBF.OBS.SOM');


%Note: all "otherpriors" and "otherpriorsunc" will be passed as follows
if isfield(CBF.OBS,'ABGB') & isempty(CBF.OBS.ABGB)==0 & CBF.OBSUNC.ABGB.unc==-9999;error('Need to prescribe CBF.AGBGunc field');end
if isfield(CBF.OBS,'SOM') & isempty(CBF.OBS.SOM)==0 & CBF.OBSUNC.SOM.unc==-9999;error('Need to prescribe CBF.SOMunc field');end

end

%other (time-invariant) constraints 
function [OPR,OPRU]=write_other_obs_constraints(CBF,OPR,OPRU);

%Mean biomass
OPR(1)=CBF.OTHER_OBS.MBiomass.mean;
OPRU(1)=CBF.OTHER_OBS.MBiomass.unc;

%Mean fire emissions
OPR(3)=CBF.OTHER_OBS.MFire.mean;
OPRU(3)=CBF.OTHER_OBS.MFire.unc;

%Mean LAI
OPR(5)=CBF.OTHER_OBS.MLAI.mean;
OPRU(5)=CBF.OTHER_OBS.MLAI.unc;

%Mean GPP
OPR(6)=CBF.OTHER_OBS.MGPP.mean;
OPRU(6)=CBF.OTHER_OBS.MGPP.unc;

%




end


