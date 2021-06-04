function ph=plotglobal(A,x,y,mapopt,lw,col)
%FUNCTION PLOTGLOBAL(A)
%function plotglobal(A,x,y,mapopt)
%
%Assumes A is a 2D array spanning the globe on eurocentric cartesian
%coordinates.
%
%Does all the crap stuff for you

%Setting the coordinates
defval('lw',1)
defval('x',[])
defval('mapopt',1)
defval('col',[0,0,0])

if isempty(x)
xres=360/size(A,2);
yres=180/size(A,1);

xvec=-180+xres/2:xres:180-xres/2;
yvec=-90+yres/2:yres:90-yres/2;
elseif length(x)<numel(x)
    
    xvec=x(1,:);
    
    yvec=y(:,1);
    
else
    xvec=x;
    yvec=y;
end


ph=imagesc(xvec,yvec,A);
set(gca,'YDir','Normal')

hold on
plotworldmap(mapopt,[],col,lw)

