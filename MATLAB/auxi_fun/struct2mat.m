function M=struct2mat(S,fldnm)
%Struct S must be solely comprised of double values.

    
    

fn=fieldnames(S);
M=[];
if nargin<2;fldnm='';end


for f=1:numel(fn);
    
    if isstruct(S.(fn{f}));
        M=[M;struct2mat(S.(fn{f}),fldnm)];
    elseif isempty(fldnm) | strcmp(fldnm,fn{f})
    M=[M;S.(fn{f})(:)];
    end
end





end