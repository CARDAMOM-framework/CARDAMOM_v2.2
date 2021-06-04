function plotworldmap(opt,rivers,col,lw)

defval('lw',1)
defval('col',[0,0,0])
%plots world map with underlying world image as main option
%option 1: PLOTS BLANK MAP
%option 2: PLOTS IMAGE SAT.

defval( 'opt',1)
defval('rivers',0)
switch opt
case 2
    %Loading image SATMAP
    %
    %Map obtained from: http://www.ham-radio-deluxe.com/Downloads/BetaKits/tabid/99/Default.aspx
    satmap=imread('phdmatlab/imag_dat/earth-living-1600-800.jpg');
    %Creating Coordinate System
    %origin = [-180,90]
    xcoords=([1:1600] - 1)/1600*360-180;
    ycoords=([800:-1:1] - 1)/800*180-90;
    imagesc(xcoords,ycoords,satmap);
    set(gca,'YDir','Normal')
    %Done with map
    %
    %Now using Matlab colorscale to overplot data
    %
    case 1
        %plotting only continental outlines
        load('CARDAMOM/DATA/MAPS/worldmap.mat');
        plot(lonlat(:,1),lonlat(:,2),'Color',col)
    case 3 
        %omit antarctica;
                load('/CARDAMOM/DATA/MAPS/worldmap.mat')
 
        lonlata=lonlat;
        plot(lonlat(:,1),lonlat(:,2),'Color',col)
                        lonlata(lonlat(:,2)>-60,:)=NaN;
        plot(lonlata(:,1),lonlata(:,2),'Color',[1,1,1]*0.9)
end

hold on


%plotting rivers
if rivers
    load riverdata
    %Plotting: Amazon(111,112,118), Congo(113,115), Ganges(88,89),
    %Niger(99,103), Nile(77,78,105, 107),Tigris(25), Euphrates(26),
    %Danube(54), Volga(50) Mississippi +Missouri (60,66,71,83), Indus(86), 
    %Ganges+brahma(87,88,89), Orinoco(104), Yangtze(84), Xi(90), Yellow(69,73)
    rivindex=[111,112,118,113,115,88,89,99,103,77,78,105,107,75,76,54,...
        50,60,66,71,69,73,83,84,86,87:89,104,90];
    for n=1:length(rivindex)
        plot(rivers(rivindex(n)).Lon,rivers(rivindex(n)).Lat,'Color',col,'Linewidth',lw)

    end
end

    
    