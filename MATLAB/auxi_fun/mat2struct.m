function [S,k]=mat2struct(M,S0,k)
%M is a vector/array of double-precision values
%S0 is a structure with pre-defined double values 

S=S0;
fn=fieldnames(S);
if nargin==2;k=0;end
for f=1:numel(fn)

    if isstruct(S.(fn{f}));
        [S2,k]=mat2struct(M,S.(fn{f}),k);
        S.(fn{f})=S2;
    else
        N=numel(S.(fn{f}));
            
    S.(fn{f})(1:N)=M(k+[1:N]);
        k=k+N;
    end
end



end