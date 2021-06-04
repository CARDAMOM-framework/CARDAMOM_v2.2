function TSA=monthly2seasonal(TS,nm,meanmedian)
%Should be rows for 2D mat
%If 3D mat, assumes 3rd dim is "seasonalized"

TS=squeeze(TS);

if total(size(TS))>numel(TS) & size(TS,1)>size(TS,2) & ndims(TS)==2;TS=TS';end



defval('nm',1) %nm = nanmean
defval('meanmedian',1);


switch ndims(TS)
    
    case 2

if meanmedian==1
       
if nm==0;
for m=1:12;TSA(:,m)=mean(TS(:,m:12:end),2);end
else
    for m=1:12;TSA(:,m)=nanmean(TS(:,m:12:end),2);end
end

elseif meanmedian==2
    
    if nm==0;
for m=1:12;TSA(:,m)=median(TS(:,m:12:end),2);end
else
    for m=1:12;TSA(:,m)=nanmedian(TS(:,m:12:end),2);end
    end

end

    case 3

if meanmedian==1
       
if nm==0;
for m=1:12;TSA(:,:,m)=mean(TS(:,:,m:12:end),3);end
else
    for m=1:12;TSA(:,:,m)=nanmean(TS(:,:,m:12:end),3);end
end

elseif meanmedian==2
    
    if nm==0;
for m=1:12;TSA(:,:,m)=median(TS(:,:,m:12:end),3);end
else
    for m=1:12;TSA(:,:,m)=nanmedian(TS(:,:,m:12:end),3);end
    end

end

end





end