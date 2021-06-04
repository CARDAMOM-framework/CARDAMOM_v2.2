function [h,x,y]=plotpdf(values,probint,minb,maxb,nobins,col,weights)
% [h,x,y]=plotpdf(values,probint,minb,maxb,nobins,col,weights)
%this function uses hist_factor to plot values as a probabo;otu density
%function!
%write more info!

%TEST: [h,x,y]=plotpdf(randn(1,ceil(rand(1)*10000)),1,-10,10,80)
%this function should always give a constant result!

if nargin==2;
   col=probint;
   clear probint
end
defval('probint',1)
%calculate order of magnitude
%oom=10.^ceil(log10(maxi(values)/mini(values)));
defval('minb',floor(mini(values)))
defval('maxb',ceil(maxi(values)))
defval('nobins',100);
defval('col',rand(1,3)*0.9)
defval('weights',1)

if minb>maxb;warning('Swapping min and max values');mx=minb;minb=maxb;maxb=mx;end

res=(maxb-minb)/nobins;
% probint/res./numel(values)

[x,y]=hist_factor(values,weights,minb,maxb,nobins);
%total of y = numel;

%normalising to give likelihood per probint
%i.e. total of y = probint/res
y=y.*probint/(res*total(y));

h=plot(x,y,'Color',col,'LineWidth',2);
