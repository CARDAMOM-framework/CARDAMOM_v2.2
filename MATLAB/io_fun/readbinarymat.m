function V=readbinarymat(file,dims,doubleint)
%V=readbinarymat(FILE,[r,c])
%
%reads binary file FILE
%FILE contains N elements to be ordered with r rows and c columns, where
%N=r x c; 
%
%function written for double precision values (e.g. MHMCMC output)
%
%if r=0, then it is set to r=N/c
%if c=0, then it is set to c=N/r

fd=fopen(file);

if nargin<3;doubleint=1;end



if doubleint==1
VV=fread(fd,inf,'real*8');
elseif doubleint==2
    VV=fread(fd,inf,'int');
end
    
%If one dimension is omitted (i.e. dims(n)==0) then it is completed based
%on the remaining requested dims and number of elements in file
if nargin<2 | isempty(dims);
    dims=[numel(VV),1];
end
for n=1:numel(dims);if dims(n)==0;dims(n)=1;dims(n)=numel(VV)/prod(dims);end;end


V=zeros(dims([2,1,3:end]));
%swapping dimensions for auto-completion
%permute(V,[2,1,3:ndims(V)]);

if numel(V)~=numel(VV);
    disp('Warning: DIMS do not match the file size');
    disp(['number of elements in binary file = ',num2str(numel(VV))]);
    disp(['number of elements in requested array  = ',num2str(numel(V))]);
    disp(['Dimensions requested = ',num2str(size(V))]);
end




V(1:end)=VV;
s=size(V);
%swapping dimensions for row, column,page order
permute(V,[2,1,3:ndims(V)]);

fclose(fd);



end