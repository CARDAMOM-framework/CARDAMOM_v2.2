function V=readbinarymat(filename,variable,name)
disp(sprintf('Saving %s ',filename));
S.(name) = variable;
save(filename,'-struct','S');
end