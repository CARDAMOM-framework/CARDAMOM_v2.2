function b=normvec(a,md);
defval('md',1)
%normalizes vector by dividing quantity by the mean of the quantity
%md 3: normalize by removing mean
%
%if a is array, then normalize each row
if length(a)<numel(a)
    for n=1:size(a,1)
        b(n,:)=normvec(a(n,:),md);
    end
else
            

    if md==1
    b=a/mean(a);
    elseif md==2;
        %coefficient of variation
        b=(a-mean(a))/std(a);
    elseif md==3;
        b=a-mean(a);
    elseif md==4;
        b=a/total(a);
    end

end
    


end