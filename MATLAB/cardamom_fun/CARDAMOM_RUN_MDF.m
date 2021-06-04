function D=CARDAMOM_RUN_MDF(CBF,MCO,cbrfile,compile,command_only,Cpath,STARTPARS)
%This function sets up a CARDAMOM Metropolis-Hastings Markov Chain Monte Carlo (MHMCMC) chain:
%INPUTS
% - CBF: a CBF structure OR a CBF filename
% - MCO: a structure containing all MHMCMC options  (default "MCO strucure
% included below)
% -"compile" (1 or 0): ensures latest version of C code is compiled
% (default = 1);
%
%Default MCO structure
%   MCO.niterations=10000;
%   MCO.printrate=1000;
%   MCO.samplerate=10;
% MCO.minstepsize=1e-4
%
%See also: CARDAMOM_WRITE_BINARY_FILEFORMAT.m
%Last modified by A.A. Bloom 2020/01/25
if nargin<2; MCO=[];end

if isnumeric(MCO) & numel(MCO)==3
    %Can also provide numerical inputs
    MCOnum=MCO;clear MCO;MCO.niterations=MCOnum(1);MCO.printrate=MCOnum(2);MCO.samplerate=MCOnum(3);
elseif isstruct(MCO) | isempty(MCO);
    
    MCOin=MCO;

if isfield(MCO,'niterations')==0;MCO.niterations=100000;end
if isfield(MCO,'printrate')==0;MCO.printrate=1000;end
if isfield(MCO,'samplerate')==0;MCO.samplerate=MCO.niterations/2000;end
if isfield(MCO,'minstepsize')==0;MCO.minstepsize=1e-5;end
if isfield(MCO,'mcmcid')==0;MCO.mcmcid=119;end
if isfield(MCO,'nadapt')==0;MCO.nadapt=100;end
disp(MCO);
if numel(fields(MCO))>6; 
    disp('*****WARNING*****')
    disp('Unused fieldnames provided in MCMC options structure "MCO"; while these won''t interfere with CARDAMOM_RUN_MDF, check for typos just in case');
    disp('******************')
    keyboard;
end

end


cbrtemp=0;cbftemp=0;

if nargin<3; a=rng; rng(mod(now*10000,2^32));cbrfile=['DUMPFILES/CBR',char(ceil(rand(1,20)*25)+96),'temp.cbr'];cbrtemp=1;rng(a);end
if nargin<6;Cpath=getenv('CARDAMOM_C_PATH');end
if nargin<5 | isempty(command_only); command_only=0;end
if nargin<4 | isempty(compile); compile=1;end
%CARDAMOM CBF file
if isstr(CBF)==0;
    cbffile=['DUMPFILES/CBF',char(ceil(rand(1,20)*25)+96),'temp.cbf'];
    CARDAMOM_WRITE_BINARY_FILEFORMAT(CBF,cbffile);
    cbftemp=1;
else
    cbffile=CBF;
end




if compile==1
CARDAMOM_COMPILE(Cpath)
end




%compilation command
%cp=unix(sprintf('gcc %s/projects/CARDAMOM_MDF/CARDAMOM_MDF.c -o %s/projects/CARDAMOM_MDF/CARDAMOM_MDF.exe -lm',Cpath,Cpath));
%ensures code is compiled





%Adding 4th option 
%INPUTS
%CARDAMOM data structure OR cardamom data file


%Outputs
%output = data


                    %make start file 
                        %make_start_file
            
            %mdf command
                       
                            command_cdea=sprintf('%s/projects/CARDAMOM_MDF/CARDAMOM_MDF.exe %s %s %d %d %d %f %d %d',Cpath,cbffile, cbrfile,MCO.niterations,MCO.printrate,MCO.samplerate,MCO.minstepsize,MCO.mcmcid,MCO.nadapt);
                       
                       %run command
                       disp('C executable command');
                       disp(command_cdea);
                       disp('******************');
if command_only==0
    
    uxid=unix(command_cdea);

                       if uxid>0; error(sprintf('%s\n ERROR!!! the unix execution of the above command failed.',command_cdea));end
%run results
D=CARDAMOM_RUN_MODEL(cbffile,cbrfile);
    D.run_mode='inverse';        


%delete files
if cbftemp==1;delete(cbffile);end
if cbrtemp==1;delete(cbrfile);end


else 
    D=command_cdea;

end




end




function [PARSALL,ANCILLARY]=read_cbr_file(filename,INFO)

%Read CBR file
PARSALL=[];

%for wildcard option; 
if isstr(filename);filename=auxifun_fullpathdir(filename);end

for n=1:numel(filename)

disp(sprintf('CHECK: .cbr file "%s" successfully read into matlab.',filename{n}));


    %number of parameter sets 
    fd=fopen(filename{n});if fd==-1;disp('WARNING: file not opened, expect error!!');end;av=fread(fd,inf,'real*8');fclose(fd);N=numel(av)/INFO.nopars;
    PARS=readbinarymat(filename{n},[N,INFO.nopars])';


if INFO.latterhalf==1;PARS=PARS(ceil(end/2)+1:end,:);end

PARSALL=[PARSALL;PARS];

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



end




