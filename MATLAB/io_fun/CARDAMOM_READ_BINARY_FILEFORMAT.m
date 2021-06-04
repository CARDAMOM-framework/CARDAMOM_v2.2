function [CBF,ANCILLARY]=CARDAMOM_READ_BINARY_FILEFORMAT(filename,INFO)
%INPUTS: "filename" and INFO (for .cbr files only)
%For .cbf format, no "INFO" structre is required.
%OUTPUTS: CBF structure, which contains all the contents of the .cbf format
%filename
%
%.cbf = (c)ARDAMOM (b)INARY (f)ILE
%.cbr = (c)ARDAMOM (b)INARY (r)ESULTS
%See also:
%"CARDAMOM_WRITE_BINARY_FILEFORMAT.m"
%s
%Last modified by A.A. Bloom 2020/05/27
%git test 7 on Feb 26 2020

if nargin==0
    error('Need to provide .cbf or .cbr file');
elseif nargin==1 
    if filename(end)=='r'; %i.e. reading cbr file
    error('Need to provide INFO struct for cbr file');
    else
        INFO=[];
    end
end


if isfield(INFO,'latterhalf')==0;INFO.latterhalf=1;end



%OPTION A: filename is CBF
if iscell(filename) | strcmp(filename(end-2:end),'cbr') | strcmp(filename(end-7:end),'cbrSTART');[CBF,ANCILLARY]=read_cbr_file(filename,INFO);
elseif strcmp(filename(end-2:end),'cbf');[CBF,ANCILLARY]=read_cbf_file(filename);end
%OPTION B: filename is CBR

end

function [CBF,BD]=read_cbf_file(filename)

%Read CBF file


%opening, reading and closing DALEC binary input file
fd=fopen(filename);
BD=fread(fd,inf,'real*8');
fclose(fd);
%can also add option to display binary data contents in friendly manner...
%OK, doing it now;
k=0;
%'Static data (100 places)')
SD=BD(k+1:k+100);
k=k+100;
%'Priors (50 places)')
PR=BD(k+1:k+50);
k=k+50;
%'Priorunc (50 places)')
PRU=BD(k+1:k+50);
k=k+50;

%'O. Priors (50 places)')
OPR=BD(k+1:k+50);
k=k+50;
%'O. Priorunc (50 places)')
OPRU=BD(k+1:k+50);
k=k+50;

%Prior information
CBF.PARPRIORS=PR;
CBF.PARPRIORUNC=PRU;
%Other constraints 
%List explicitly here
CBF=read_other_obs_constraints(CBF,OPR,OPRU);


%ID (not used)
CBF.ID=SD(1);
%Latitude
CBF.LAT=SD(2);
%Number of days
CBF.nodays=SD(3);
%number of met fields
CBF.nomet=SD(4);
%Number of observation streams
CBF.noobs=SD(5);
%Implement EDCs
CBF.EDC=SD(6);
%EDC diagnostic
CBF.EDCDIAG=SD(7);



%MCMC start searching EDCs from anywhere (1) or from prescribed starting
%point(0). this is outdated - consider deleting
CBF.rc_random_search=SD(11)==1;BD(11)=-9999;

%NEE IAV options
CBF=read_obs_uncertainty_fields(CBF,SD,OPRU);


TEMPDATA=zeros(CBF.nomet+CBF.noobs,CBF.nodays);
TEMPDATA(1:end)=BD(k+1:k+numel(TEMPDATA));

%All met data
CBF.MET=TEMPDATA(1:CBF.nomet,1:CBF.nodays)';
%All observations (if any) 
%1st column = GPP
%2nd column = LAI
%3rd column = NEE
CBFOBS=TEMPDATA(CBF.nomet+1:end,1:CBF.nodays)';
CBF=define_cbf_obs_fields(CBF,CBFOBS);

disp(sprintf('CHECK: .cbf file "%s" successfully read into matlab.',filename));

%Removing redundant fields
CBF=rmfield(CBF,'noobs');
% CBF=rmfield(CBF,'nomet');
% CBF=rmfield(CBF,'nodays');

%Retaining "OTHERPRIORS" for now
CBF.RAW.OTHERPRIORS=OPR;
CBF.RAW.OTHERPRIORSUNC=OPRU;
CBF.RAW.info='Raw inputs/outputs as stored in CBF binary structure';
CBF.RAW.details='For completeness & development purpose only; When re-writing CBF to file, these are over-written by CBF.OBS, etc.';



end


%Time-resolved observations
function CBF=define_cbf_obs_fields(CBF,CBFOBS);


obsnames={'GPP','LAI','NBE','ABGB','ET','EWT','BAND1','BAND2','BAND3','BAND4','SOM','CH4'};
for n=1:numel(obsnames);%CBF.noobs;
       if size(CBFOBS,2)>=n && sum(CBFOBS(:,n)~=-9999);
            CBF.OBS.(obsnames{n})=CBFOBS(:,n);
       else
                       CBF.OBS.(obsnames{n})=[];
       end
end
    

CBF.OBSinfo.LAI='LAI data requirements: should be +ve definite (>0); assumed uncertainty for LAI timeseries is factor of 2; future versions will include user-defined uncertainty options';
CBF.OBSinfo.uncertainty_factors='Uncertainty structures for positive-definite quantities (GPP, LAI, ET), are prescribed as uncertainty factors (by default); uncertainty factors should be > 1. \n For example: a 1-sigma range for 100 uncertainty factor 2 = 100/2 - 100*2 = 50 - 200 ';
end
%Uncertainties for time-resolved observations
function CBF=read_obs_uncertainty_fields(CBF,SD,OPRU);

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
CBF.OBSUNC.NBE.annual_unc=SD(14);
CBF.OBSUNC.NBE.seasonal_unc=SD(16);
CBF.OBSUNC.NBE.info='Single point (default = 0.5, must be >0) and annual (annual_unc) NBE uncertainty [gC/m2/d]';

%ET 
CBF.OBSUNC.ET.annual_unc=SD(15);
CBF.OBSUNC.ET.unc=SD(17);
CBF.OBSUNC.ET.obs_unc_threshold=SD(22);%Default = 0.1
CBF.OBSUNC.ET.info=sprintf('Single point uncertainty factor (default = */ +2, must be >1) and annual (annual_unc) ET uncertainty [factor]. \nDefault obs threshold = 0.1 mm/day: this ensures log-tranformed model ET values are insensitive to ET<0.1');

%EWT
CBF.OBSUNC.EWT.annual_unc=SD(18);
CBF.OBSUNC.EWT.unc=SD(19);
CBF.OBSUNC.EWT.info=sprintf('Single point (default = 50) and annual (N/A yet) EWT uncertainty [mm]');

%GPP 
CBF.OBSUNC.GPP.annual_unc=SD(20);
CBF.OBSUNC.GPP.unc=SD(21);
CBF.OBSUNC.GPP.info=sprintf('Single point uncertaint factor (default = */ +2, must be >1) and annual (annual_unc) GPP uncertainty [gC/m2/d]\n default obs unc threshold is 0.1gC/m2/d: this ensures log-tranformed model GPP values are insensitive to GPP<0.1');
%gpp abs: option on treating GPP constraint as absolute or relative
CBF.OBSUNC.GPP.gppabs=SD(8);
CBF.OBSUNC.GPP.obs_unc_threshold=SD(23);%Default = 0.1
CBF.OBSUNC.GPP.gppabs_info='Set to "1" for GPP data, set to "0" for SIF data"';


%CH4
CBF.OBSUNC.CH4.annual_unc=SD(25);
CBF.OBSUNC.CH4.unc=SD(26);
CBF.OBSUNC.CH4.info=sprintf('Single point uncertaint factor (default = */ +2, must be >1) and annual (annual_unc) CH4 uncertainty [mgC/m2/d]\n default obs unc threshold is 0.1mgC/m2/d: this ensures log-tranformed model CH4 values are insensitive to CH4<0.1');
CBF.OBSUNC.CH4.obs_unc_threshold=SD(27);%Default = 0.01


%Time-resolved SOM uncertainty
CBF.OBSUNC.SOM.unc=OPRU(8);
CBF.OBSUNC.SOM.info=sprintf('Single point uncertainty on time-resolved SOM, CBF.OBS.SOM');

%Time-resolved biomass uncertaiknty
CBF.OBSUNC.ABGB.unc=OPRU(2);
CBF.OBSUNC.ABGB.info=sprintf('Uncertainty *factor* on time-resolved biomass CBF.OBS.ABGB');



end
%other (time-invariant) constraints 
function CBF=read_other_obs_constraints(CBF,OPR,OPRU)
%Mean biomass
CBF.OTHER_OBS.MBiomass.mean=OPR(1);
CBF.OTHER_OBS.MBiomass.unc=OPRU(1);
CBF.OTHER_OBS.MBiomass.info=sprintf('Mean total (above-& below- ground) biomass at time t=0\n Uncertainty distribution is log-normal distribution');

%Mean fire emissions
CBF.OTHER_OBS.MFire.mean=OPR(3);
CBF.OTHER_OBS.MFire.unc=OPRU(3);
CBF.OTHER_OBS.MFire.info=sprintf('Mean total (above-& below- ground) fire emissions for whole simulation \n Use +ve value for uncertainty factor  (log-gaussian distribution) & -ve value for absolute uncertainty (naussian distribution)');
%Mean LAI
CBF.OTHER_OBS.MLAI.mean=OPR(5);
CBF.OTHER_OBS.MLAI.unc=OPRU(5);
CBF.OTHER_OBS.MLAI.info=sprintf('Mean total LAI for whole simulation \n Use +ve value for uncertainty factor  (log-gaussian distribution) & -ve value for absolute uncertainty (gaussian distribution)\n If any CBF.OBS.LAI values are presribed, then mean LAI will be used to constrain these timesteps only.'  );
%Mean GPP
CBF.OTHER_OBS.MGPP.mean=OPR(6);
CBF.OTHER_OBS.MGPP.unc=OPRU(6);
CBF.OTHER_OBS.MGPP.info=sprintf('Mean total (above-& below- ground) GPP for whole simulation \n Use +ve value for uncertainty factor  (log-gaussian distribution) & -ve value for absolute uncertainty (naussian distribution)');




end



function [PARSALL,ANCILLARY]=read_cbr_file(filename,INFO)

%Read CBR file
PARSALL=[];
CHAINALL=[];

%for wildcard option; 
if isstr(filename);filename=auxifun_fullpathdir(filename);end

for n=1:numel(filename)

disp(sprintf('CHECK: .cbr file "%s" successfully read into matlab.',filename{n}));


if isempty(dir(filename{n}))==0
    %number of parameter sets 
    fd=fopen(filename{n});if fd==-1;disp('WARNING: file not opened, expect error!!');end;av=fread(fd,inf,'real*8');fclose(fd);N=numel(av)/INFO.nopars;
    PARS=readbinarymat(filename{n},[N,INFO.nopars])';


if INFO.latterhalf==1;PARS=PARS(ceil(end/2)+1:end,:);end
else
    PARS=[];
end

PARSALL=[PARSALL;PARS];
CHAINALL=[CHAINALL;n*ones(size(PARS,1),1)];

end

% %reading step files in case there
% for f=1:numel(filename)
%     ancfile=[filename{f},'STEP'];
%     if isempty(dir(ancfile))==0;
%         fopen(ancfile);av=fread(fd,inf,'real*8');fclose(fd);N=numel(av)/INFO.nopars;
%             STEP{f}=readbinarymat(ancfile,[N,INFO.nopars])';
%     else
%         STEP{f}=[];
%     end
% end

%reading start files in case there
for f=1:numel(filename)
    ancfile=[filename{f},'START'];
    if isempty(dir(ancfile))==0;
        fopen(ancfile);av=fread(fd,inf,'real*8');fclose(fd);N=numel(av)/INFO.nopars;
            START{f}=readbinarymat(ancfile,[N,INFO.nopars])';
    else
        START{f}=[];
    end
end


ANCILLARY.START=START;
ANCILLARY.STEP=[];
ANCILLARY.chainid=CHAINALL;





end


%***************

