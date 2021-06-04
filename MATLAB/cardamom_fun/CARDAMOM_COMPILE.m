function CARDAMOM_COMPILE(Cpath,opt)



if nargin==0 | isempty(Cpath);
    Cpath=getenv('CARDAMOM_C_PATH');
end

if nargin<2
        opt=1;
end

   disp(sprintf('CARDAMOM compilation based on scripts in %s',Cpath));
compile_command_run_model=sprintf('gcc %s/projects/CARDAMOM_GENERAL/CARDAMOM_RUN_MODEL.c -o %s/projects/CARDAMOM_GENERAL/CARDAMOM_RUN_MODEL.exe -lm',Cpath,Cpath);
compile_command_run_mdf=sprintf('gcc %s/projects/CARDAMOM_MDF/CARDAMOM_MDF.c -o %s/projects/CARDAMOM_MDF/CARDAMOM_MDF.exe -lm',Cpath,Cpath);
   compile_command_run_assemble=sprintf('gcc %s/projects/CARDAMOM_GENERAL/CARDAMOM_ASSEMBLE_MODELS.c -o %s/projects/CARDAMOM_GENERAL/CARDAMOM_ASSEMBLE_MODELS.exe -lm',Cpath,Cpath);

if opt>=1;
cp(1)=unix(compile_command_run_model);
cp(2)=unix(compile_command_run_mdf);
cp(3)=unix(compile_command_run_assemble);
if opt==2
    %Debug compilation as well
    %gcc /Users/abloom/EDI/C/projects/CARDAMOM_MDF/CARDAMOM_MDF.c -g -o /Users/abloom/EDI/C/projects/CARDAMOM_MDF/CARDAMOM_MDF.exe -lm
compile_command_run_model_debug=sprintf('gcc %s/projects/CARDAMOM_GENERAL/CARDAMOM_RUN_MODEL.c -g -o %s/projects/CARDAMOM_GENERAL/CARDAMOM_RUN_MODEL_debug.exe -lm',Cpath,Cpath);
compile_command_run_mdf_debug=sprintf('gcc %s/projects/CARDAMOM_MDF/CARDAMOM_MDF.c -g -ggdb3 -o %s/projects/CARDAMOM_MDF/CARDAMOM_MDF_debug.exe -lm',Cpath,Cpath);
    cpd(1)=unix(compile_command_run_model_debug);
cpd(2)=unix(compile_command_run_mdf_debug);
end

if sum(cp)>0;
    warning('Step 1: Warning - CARDAMOM_RUN_MODEL.exe or  CARDAMOM_RUN_MDF.exe did not compile properly');
else

        disp('CARDAMOM code successfully compiled');
end


elseif opt==0;
    disp('Compilation command:')
    disp(compile_command_run_assemble)
    disp(compile_command_run_model)
    disp(compile_command_run_mdf);
end

%Final step: run "assemble" command to produce text files.

%unix(sprintf('./%s/projects/CARDAMOM_GENERAL/CARDAMOM_ASSEMBLE_MODELS.exe',Cpath));




end
