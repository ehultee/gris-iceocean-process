clear; close all;

% load ismip bed
bedfile = 'ismip_bed.nc';
x = double(ncread(bedfile,'x'));
y = double(ncread(bedfile,'y'));
bed = double(ncread(bedfile,'bed'));

% read convex hull (just for plotting)
T = readtable('convex_hull.csv');

% working plot
figure(); hold on;
pcolor2(x,y,bed');
axis xy equal;
plot(T.xb,T.yb,'k','linewidth',1);
colorbar; clim([-1000,0]);

% deep ocean pixels (what we're looking for connection to)
do = double(bed<-2000);

% get effective depths
z_eff = NaN*bed;
z_eff_res = 50;
z_eff_min = -ceil(abs(min(bed(:)))/z_eff_res)*z_eff_res;
zlevels = 0:-z_eff_res:z_eff_min;
for k=1:length(zlevels)
    bw = double(bed<=zlevels(k));
    regions = bwconncomp(bw);
    for j=1:regions.NumObjects
        inds = regions.PixelIdxList{j};
        if any(do(inds))
            z_eff(inds) = zlevels(k);
        end
    end
end

% give above sea level points an effective depth of 100 m
z_eff(bed>0) = 100;

% assign a minimum effective depth of -2000 m
% this prevents problems when the cmip model does not have sufficiently
% deep points, and effective depths deeper than -2000 m are very
% unlikely to be relevant to the ice sheet simulations
z_eff(z_eff<-2000) = -2000;

% plot effective depths
figure(); hold on;
pcolor2(x,y,z_eff');
axis xy equal;
plot(T.xb,T.yb,'k','linewidth',1);
colorbar; clim([-1000,0]);

% create netcdf
output_file = 'z_eff.nc';
if exist(output_file, 'file'), delete(output_file); end

% initialise variables in netcdf
nccreate(output_file,'x','Dimensions',{'x',length(x)},'Datatype','single');
nccreate(output_file,'y','Dimensions',{'y',length(y)},'Datatype','single');
nccreate(output_file,'z_eff','Dimensions',{'x',length(x),'y',length(y)},'Datatype','single','deflatelevel',9);

% write variables
ncwrite(output_file,'x',single(x));
ncwrite(output_file,'y',single(y));
ncwrite(output_file,'z_eff',single(z_eff));