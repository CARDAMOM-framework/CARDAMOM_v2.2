function MA=CARDAMOM_MODEL_LIBRARY(ID,MA,reread)
%MA=CARDAMOM_MODEL_LIBRARY(ID)
%
%INPUTS: CARDAMOM model ID (integer) OR structure with ID field. 
%This function reads C code and extracts model attributes
%[Cpath]/projects/CARDAMOM_GENERAL/CARDAMOM_MODEL_LIBRARY.c
%OUTPUTS: "MA" structure with CARDAMOM model attributes:
%Specifically, the structure contains
% MA.nopools: number of model pools
%MA.nopars: number of model parameters
%MA.nofluxes: number of stored model fluxes
%
%Example: 
%   MA=CARDAMOM_MODEL_LIBRARY(2);
%
%Last modified by A.A. Bloom 2020 June 17
Cpath=getenv('CARDAMOM_C_PATH');

%C script
if nargin<2;MA=[];end
if nargin<3;reread=0;end

%if isstruct(ID);MA=ID;ID=MA.ID;end
if isstruct(ID);ID=ID.ID;end
dumpfile=sprintf('DUMPFILES/CARDAMOM_MODEL_LIBRARY_ID=%i.mat',ID);
%TO DO: access new model types here (805 etc.) eventually making access to
%CARDAMOM_MODEL_LIBRARY.c obsolete



if reread==1
%*****Step 0. Re-run CARDAMOM_ASSEMBLE_MODELS.exe
unix(sprintf('./%s/projects/CARDAMOM_GENERAL/CARDAMOM_ASSEMBLE_MODELS.exe %s',Cpath,Cpath));
end

%*****Step 1. find path for model info*****
filename=sprintf('%s/projects/CARDAMOM_MODELS/DALEC/DALEC_%i/MODEL_INFO_%i.c',Cpath,ID,ID);
stack_filename=sprintf('%s/projects/CARDAMOM_MODELS/DALEC/DALEC_%i/dalec_%i_pars.txt',Cpath,ID,ID);

%filename=sprintf('%s/projects/CARDAMOM_GENERAL/CARDAMOM_MODEL_LIBRARY.c',Cpath);
filedir=dir(filename);
if isempty(filedir);error('Warning: model version not currently supported by CARDAMOM');end
dumpfiledir=dir(dumpfile);

if isempty(dir(dumpfile)) | reread==1 | dumpfiledir.datenum<filedir.datenum
    disp('Re-reading file...');

D=importdata(filename,'');


%searching for model attribute" line zero"
%this line contains exactly '/*IDxMA*/' where x is the model ID number
%line0=find(strcmp(D,sprintf('/*ID%dMA*/',ID)));


for n=1:size(D,1);
    %reading strings that read as follows
    %DATA->nopools=7;
    if strncmp(D{n},'DALECmodel.no',12)
        %E.g. evaluates  'DALECmodel.nopools=7;'
    eval(D{n}(1:find(D{n}==';',1)));
    end
end



%Overwriting with stack file
if isfile(stack_filename)
D=importdata(stack_filename);
for n=1:size(D);
    eval(D{n});
end




else
%Find parameter info here *******
parfilename=sprintf('%s/projects/CARDAMOM_MODELS/DALEC/DALEC_%i/PARS_INFO_%i.c',Cpath,ID,ID);


%repeat the same with model parameter value ranges
%TO DO: need to set these up for other models too


 D=importdata(parfilename,'');

 pn=1;
 for n=1:size(D,1)
     linestr=D{n};
     if numel(linestr)>21 & strcmp(linestr(1:12),'CARDADATA->p')==1  
     linestr(linestr=='-')='';
     linestr(linestr=='>')='.';
     b1=find(linestr=='[');
     b2=find(linestr==']');
     linestr(b1)='(';
     linestr(b2)=')';
     lnum=num2str(str2num(linestr(b1+1:b2-1))+1);
     linestr=[linestr(1:b1),lnum,linestr(b2:end)];
     %Find first ";"
     eval(linestr(1:find(linestr==';',1)));
     %Storing parameter name based on line above parmin
     if strcmp(linestr(1:16),'CARDADATA.parmin');
         if strcmp(D{n-1}(1:2), '/*')
         parname{pn}=D{n-1};
         else
             parname{pn}='';
         end
         pn=pn+1;
     end
     
     end
 end
% 
% 
 MA.parmin=CARDADATA.parmin;
 MA.parmax=CARDADATA.parmax;
 MA.parname=parname;
 
end


%Store fields regardless
for f=fieldnames(DALECmodel)';
            MA.(f{1})=DALECmodel.(f{1});
end


 
 
 
 %Random sample
 MA.parrand=@(N) logrand(MA.parmin,MA.parmax,N);
 MA.par2nor=@(pars) log(pars./repmat(MA.parmin,[size(pars,1),1]))./log(repmat(MA.parmax./MA.parmin,[size(pars,1),1]));
 MA.nor2par=@(npars) repmat(MA.parmin,[size(npars,1),1]).*(repmat(MA.parmax,[size(npars,1),1])./repmat(MA.parmin,[size(npars,1),1])).^npars;

save(dumpfile,'MA');
else
    MAall=MA;
   load(dumpfile,'MA');
end

MA.ID=ID;

end



function MA=CARDAMOM_MODEL_LIBRARY_OLD(ID,MA,reread)
%MA=CARDAMOM_MODEL_LIBRARY(ID)
%
%INPUTS: CARDAMOM model ID (integer) OR structure with ID field. 
%This function reads C code and extracts model attributes
%[Cpath]/projects/CARDAMOM_GENERAL/CARDAMOM_MODEL_LIBRARY.c
%OUTPUTS: "MA" structure with CARDAMOM model attributes:
%Specifically, the structure contains
% MA.nopools: number of model pools
%MA.nopars: number of model parameters
%MA.nofluxes: number of stored model fluxes
%
%Example: 
%   MA=CARDAMOM_MODEL_LIBRARY(2);
%
%Last modified by A.A. Bloom 2018/02/11 



Cpath='C';

%C script
if nargin<2;MA=[];end
if nargin<3;reread=0;end

if isstruct(ID);MA=ID;ID=MA.ID;end
dumpfile=sprintf('DUMPFILES/CARDAMOM_MODEL_LIBRARY_ID=%i.mat',ID);
%TO DO: access new model types here (805 etc.) eventually making access to
%CARDAMOM_MODEL_LIBRARY.c obsolete
filename=sprintf('%s/projects/CARDAMOM_GENERAL/CARDAMOM_MODEL_LIBRARY.c',Cpath);
filedir=dir(filename);
dumpfiledir=dir(dumpfile);

if isempty(dir(dumpfile)) | reread==1 | dumpfiledir.datenum<filedir.datenum
    disp('Re-reading file...');

D=importdata(filename);


%searching for model attribute" line zero"
%this line contains exactly '/*IDxMA*/' where x is the model ID number
line0=find(strcmp(D,sprintf('/*ID%dMA*/',ID)));


for n=1:3
    %reading strings that read as follows
    %DATA->nopools=7;
    linestr=D{line0+n};
    linestr(linestr=='-')='';
    linestr(linestr=='>')='.';
    eval(linestr)
end

for f=fieldnames(DATA)';
            MA.(f{1})=DATA.(f{1});
end


%repeat the same with model parameter value ranges
%TO DO: need to set these up for other models too
if ID==1
filename=sprintf('%s/projects/DALEC_CODE/DALEC_CDEA/PARS_INFO_CDEA.c',Cpath);

D=importdata(filename);
for n=1:size(D,1)
    linestr=D{n};
    if numel(linestr)>21 & strcmp(linestr(1:12),'CARDADATA->p')==1  
    linestr(linestr=='-')='';
    linestr(linestr=='>')='.';
    b1=find(linestr=='[');
    b2=find(linestr==']');
    linestr(b1)='(';
    linestr(b2)=')';
    lnum=num2str(str2num(linestr(b1+1:b2-1))+1);
    linestr=[linestr(1:b1),lnum,linestr(b2:end)];
    eval(linestr)
    end
end


MA.parmin=CARDADATA.parmin;
MA.parmax=CARDADATA.parmax;

end;
save(dumpfile,'MA');
else
    MAall=MA;
   load(dumpfile,'MA');
end

MA.ID=ID;

end



