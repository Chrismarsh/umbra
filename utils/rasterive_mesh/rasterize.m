
files=dir('*.mat');

for f=1:max(size(files));
    
fname=files(f).name;
    
load(fname)

fprintf('%s\n',fname)

X=mxDomain(:,1);
Y=mxDomain(:,2);
[bbx,bby,~,~]=minboundrect(X,Y);

% hold on
% plot(bbx,bby,'color','r');



%bottom left most point
    [val,i]=min(bbx);
    xll=val;
    yll=bby(i);
%     plot(xll,yll,'o'); 
%bottom right most pt
    [val,i]=max(bbx);
	 xlr = val;
	 ylr = bby(i);
%      plot(xlr,ylr,'o'); 
%top right most pt
    [val,i]=max(bby);
	 yur = val;
	 xur = xlr;
%      plot(xur,yur,'o'); 
%top left most pt
    [val,i]=max(bby);
	 yul = val;
	 xul = xll;
%     plot(xul,yul,'o'); 

cell_size = 2;

num_rows= ceil( (yul-yll)/cell_size)
num_cols= ceil( abs(xll-xlr)/cell_size)

[mx,my]=meshgrid(xll:cell_size:xlr,yll:cell_size:yul);
Z=gridtrimesh(tri,mxDomain,mx,my,shadows);

Z=bwmorph(Z,'majority','inf');
% imwrite(flipud(Z),strcat('model_',fname(1:end-4)),'tif');

% header=[num_cols+1;num_rows+1;xll;yll;cell_size;-9999];
% SaveAsciiRaster(flipud(Z),header,strcat('model_',fname(1:end-4)));
end
% imshow(Z)
% 


