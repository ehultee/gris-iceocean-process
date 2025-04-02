% script to extend bedmachine far enough to cover ismip6 grid
clear; close all;

% load bedmachine (CRS EPSG:3413)
bmfile = '~/Documents/BedMachineGreenland-v5.nc';
bm.x = double(ncread(bmfile,'x'));
bm.y = double(ncread(bmfile,'y'));
bm.bx = [min(bm.x),max(bm.x),max(bm.x),min(bm.x),min(bm.x)];
bm.by = [min(bm.y),min(bm.y),max(bm.y),max(bm.y),min(bm.y)];
bm.bed = double(ncread(bmfile,'bed'));
bm.res = 150;

% define ismip6 grid (in EPSG:3413)
ismip.xg = -0.72e6:1e3:0.96e6;
ismip.yg = -3.45e6:1e3:-0.57e6;
ismip.xc = 0.5*(ismip.xg(1:end-1)+ismip.xg(2:end));
ismip.yc = 0.5*(ismip.yg(1:end-1)+ismip.yg(2:end));
% ismip6 box for plotting
ismip.bx = [ismip.xg(1),ismip.xg(end),ismip.xg(end),ismip.xg(1),ismip.xg(1)];
ismip.by = [ismip.yg(1),ismip.yg(1),ismip.yg(end),ismip.yg(end),ismip.yg(1)];

% since bedmachine does not cover a wide enough grid, augment with gebco
% read GEBCO wider bathymetry (in lat-lon)
% downloaded from https://download.gebco.net/
% after selected lat bounds 57 to 90 and lon bounds -100 to 15
gebcofile = '~/Documents/GEBCO/gebco_2024_n90.0_s57.0_w-100.0_e15.0.nc';
gb.lat = ncread(gebcofile,'lat');
gb.lon = ncread(gebcofile,'lon');
gb.bed = ncread(gebcofile,'elevation');
% subsample to make computation faster
% (a degree of lat is longer than a degree of lon so the differing
% subsampling makes it come out evenly spaced once in EPSG:3413)
ss_lat = 5;
ss_lon = 20;
gb.lon = gb.lon(1:ss_lon:end);
gb.lat = gb.lat(1:ss_lat:end);
gb.bed = gb.bed(1:ss_lon:end,1:ss_lat:end);
% put GEBCO into EPSG:3413
[gb.LAT,gb.LON] = meshgrid(gb.lat,gb.lon);
[gb.X,gb.Y] = projfwd(projcrs(3413),gb.LAT,gb.LON);

% % working plot
% figure(); hold on;
% scatter(gb.X(:),gb.Y(:),30,gb.bed(:),'filled');
% plot(ismip.bx,ismip.by,'k');
% plot(bm.bx,bm.by,'r');
% axis equal;

% create extended bedmachine grid
buffer = 50e3;
bme.xlims = [min(ismip.xg)-buffer,max(ismip.xg)+buffer];
bme.ylims = [min(ismip.yg)-buffer,max(ismip.yg)+buffer];
x1 = fliplr(min(bm.x)-bm.res:-bm.res:bme.xlims(1))';
x2 = (max(bm.x)+bm.res:bm.res:bme.xlims(2))';
bme.x = [x1;bm.x;x2];
y1 = fliplr(bm.y(1)+bm.res:bm.res:bme.ylims(2))';
y2 = (min(bm.y)-bm.res:-bm.res:bme.ylims(1))';
bme.y = [y1;bm.y;y2];
% create extended bed array and fill with the original bedmachine
bme.bed = NaN(length(bme.x),length(bme.y));
bme.bed([length(x1)+1:length(x1)+length(bm.x)],[length(y1)+1:length(y1)+length(bm.y)]) = bm.bed;

% restrict gebco to the points we need for interpolation
% must be inside ismip6 grid (+buffer)
inds = find(gb.X(:)>=bme.xlims(1)-buffer & gb.X(:)<=bme.xlims(2)+buffer...
            & gb.Y(:)>=bme.ylims(1)-buffer & gb.Y(:)<=bme.ylims(2)+buffer);
Xi = gb.X(inds);
Yi = gb.Y(inds);
bedi = gb.bed(inds);
% must be outside original bedmachine grid (-buffer)
inds = find(Xi>=min(bm.x)+buffer & Xi<=max(bm.x)-buffer...
            & Yi>=min(bm.y)+buffer & Yi<=max(bm.y)-buffer);
Xi(inds) = [];
Yi(inds) = [];
bedi(inds) = [];

% % working plot
% figure(); hold on;
% scatter(Xi,Yi,30,bedi,'filled');
% plot(ismip.bx,ismip.by,'k');
% plot(bm.bx,bm.by,'r');
% axis equal;

% do the interpolation (takes about 5 minutes)
f = scatteredInterpolant(Xi,Yi,double(bedi),'linear','none');
tic
for i=1:length(bme.x)
    inds = find(isnan(bme.bed(i,:)));
    if ~isempty(inds)
        bme.bed(i,inds) = f(bme.x(i)*ones(length(inds),1),bme.y(inds));
    end
end
toc

% final plot
figure(); hold on;
imagesc(bme.x,bme.y,bme.bed'); axis xy equal;
plot(ismip.bx,ismip.by,'k');
plot(bm.bx,bm.by,'r');

% create netcdf
output_file = 'bedmachine_extended.nc';
if exist(output_file, 'file'), delete(output_file); end

% initialise variables in netcdf
nccreate(output_file,'x','Dimensions',{'x',length(bme.x)},'Datatype','single');
nccreate(output_file,'y','Dimensions',{'y',length(bme.y)},'Datatype','single');
nccreate(output_file,'bed','Dimensions',{'x',length(bme.x),'y',length(bme.y)},'Datatype','single','deflatelevel',9);

% write variables
ncwrite(output_file,'x',single(bme.x));
ncwrite(output_file,'y',single(bme.y));
ncwrite(output_file,'bed',single(bme.bed));



