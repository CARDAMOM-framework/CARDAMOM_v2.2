function aperc=percentile(a,perc,pdim)
%function [aperc]=percentile(a,perc)
%
%please enter value as percentage (not fraction) e.g. 5%,10%, 95% are 5 10
%95.

defval('pdim',1)
%flipping for simplicity


if size(a,1)==1 & ndims(a)==2;a=a';end

aperc=perc*0;

%removing non-finite
a=a(isfinite(a));

la=size(a,pdim)-1;
a=sort(a,pdim);



dimstr='';
for n=1:ndims(a)
    if n==pdim
        dimstr=[dimstr,'round(la*perc/100)+1'];
    else
        dimstr=[dimstr,':'];
    end
    if n<ndims(a);
      dimstr=[dimstr,','];
    end
end



%only for one percentile for now...

eval(['aperc=a(',dimstr,');']);


