function ffps=auxifun_fullpathdir(fnamesearch)

%same as "dir", but does all the searching and returns full filepaths
a=dir(fnamesearch);

pathstr= fileparts(fnamesearch);

for n=1:numel(a);
    ffps{n}=fullfile(pathstr,a(n).name);
    disp(ffps{n});
end



end