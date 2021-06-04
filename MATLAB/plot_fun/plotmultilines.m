function [mh,ch,fh]=plotmultilines(x,ymat,mmopt,CMP,linecol)
%[mh,ch]=plotmultilines(x,ymat,mmopt,CMP,linecol)
%works just like plot, but plots mean with filled-in confidence ranges.
%uses all rows of ymat to make line y
%INPUTS:
%x - 1xN vector of x coordinates (set to y if only one input is assigned) 
%y - MxN array of y coordinate ensemble with M members (e.g. an ensemble of M lines)
%mmopt - plots a line of the ensemble mean if set to zero (default) or
%plots nth member if set to mmopt=n
%fullcol = the colour set for the 25% confidence range (default = red)
%backcol = the colour set for the 95% confidence range (default = light
%pink)
%linecol = the colour set for the ensemble mean, or equivalent (default =
%darker version of fullcol).
%
%OUTPUTS = handle to the ensemble mean line (or equivalent).
%
%NOTE: Function can be slow for very large ensembles.
%
%
%EXAMPLE: plot an ensemble of lines with random noise
%x=1:50; ymat=repmat(x*0.02,[200,1]); ymat=ymat+randn(size(ymat));
%plotmultilines(x,ymat)
%
%Original version by A. Anthony Bloom on 22 Oct 2013
%University of Edinburgh
%a.bloom@ed.ac.uk



%exports handle for mean (or user defined) line

defval('mmopt',0);
defval('fullcol',[1,0,0])
defval('CMP',[1,1,1;    0.9647    0.5686    0.0392])
defval('linecol',[])
defval('ymat',[])

%making x ymat if nothing in ymat
if isempty(ymat)==1;
    ymat=x;
    x=1:size(ymat,2);
    
end
nc=size(ymat,2);


%establishing colors, etc. here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of confidence ranges
CRN=20;

%confidence ranges
mpr=0+100/(CRN*2):100/CRN:100;
    
%colors for each confidence  range
COLORS=flipud(interp1(0:100/(size(CMP,1)-1):100,CMP,mpr));
% %whitening correction
% CORRECTION=((1-COLORS).^2).*[mpr',mpr',mpr']/250;
% %graying correction
% CORRECTION=(repmat(mean(COLORS,2),[1,3])-COLORS).*repmat(mpr'/500,[1,3]);

%COLORS=COLORS;
%(COLORS+CORRECTION).^repmat(1-mpr'/100,[1,3]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    

    yprc=zeros(2,nc,length(mpr));
        for nc=1:size(ymat,2)
        yprc(1,nc,:)=percentile(ymat(:,nc)',[50-mpr/2])';
        yprc(2,nc,:)=percentile(ymat(:,nc)',[50+mpr/2])';
        end
        
        hold on
    
    
  for nn=length(mpr):-1:1 

    %plotting percentile polygons
    fh(nn)=fill([x,x(end:-1:1)],[yprc(1,:,nn),fliplr(yprc(2,:,nn))],COLORS(nn,:));
    set(fh(nn),'Linestyle','none')
 end


 
 %plotting either mean or user-defined line on top of precentile polygons.
 if mmopt==0 & isempty(linecol)==0
     %mean
     mh=plot(x,median(ymat,1),'Color',linecol,'Linewidth',2);
 elseif mmopt>0 & isempty(linecol)==0
 %plot nth curve as main curve
     mh=plot(x,ymat(mmopt,:),'Color',linecol)   ;    
 else
     mh=[];
 end
 
 
 %COLORBAR
 colormap(flipud(COLORS));
 ch=colorbar;
  perclow=fliplr((100-mpr)/2);
  perchigh=100-perclow;
  caxis([min(perclow),max(perclow)]+diff(perclow(1:2))*[0,1])

  %subsampling perclow
  perclow=5:10:45;perchigh=100-[5:10:45];
  
  set(ch,'Ytick',perclow)
  yl=ylabel(ch,'confidence range');set(yl,'FontSize',get(gca,'FontSize'));
  for n=1:length(perclow)
      cblabs{n}=[num2str(perchigh(n)-perclow(n),'%2.1f'),'%'];
  end
    set(ch,'ytickLabel',cblabs)
    set(ch,'ydir','reverse')
    set(ch,'FontSize',get(gca,'FontSize'));
  %OK done with colorbar

end


