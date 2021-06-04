function [bins,pdf]=hist_factor(values,wfactor,pdfmin,pdfmax,nobins)
%function [bins,pdf]=hist_factor(values,factor,pdfmin,pdfmax,nobins)
%This is a weighted histogram function
%FACTOR determines what weight is allocated to each element in VALUES
%
%NOTE for likelihood at interval OTHER than the interval prescribed, then
%please multiply accordingly (e.g. PDF at Interval=B = PDF at
%Interval=A*B/A


defval('pdfmin',min(min(values)))
defval('wfactor',values(1).^0)
defval('pdfmax',max(max(values)))
defval('nobins',10)

pdfrange=pdfmax-pdfmin;
%re-scaling all values from 0 to nobins

V=(values-pdfmin) * nobins /(pdfrange);
Vbins=0.5:1:nobins-0.5;
%rounding to middle intervals
V=ceil(V);
pdf=Vbins*0;

if numel(wfactor)==1;
    pdf=hist(V(1:end),1:nobins)*wfactor(1);
elseif numel(wfactor)>1;
    for n=1:nobins
        pdf(n)=sum(wfactor(V==n));
    end
end
bins=Vbins*(pdfrange)/nobins+pdfmin;


end
