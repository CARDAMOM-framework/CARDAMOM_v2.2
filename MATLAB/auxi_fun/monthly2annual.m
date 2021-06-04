function TSA=monthly2annual(TS,dim)
%Should be rows for 2D mat

%default = 2 dims
defval('dim',ndims(TS))

TSA=0;

if dim==2;
TS=squeeze(TS);
if total(size(TS))>numel(TS) & size(TS,1)>size(TS,2);TS=TS';end

for m=1:12;TSA=TSA+TS(:,m:12:end)/12;end

elseif dim==3
for m=1:12;TSA=TSA+TS(:,:,m:12:end,:,:)/12;end

else 
    error('number of dims not supported, edit code accordingly!')
end
    
    
end