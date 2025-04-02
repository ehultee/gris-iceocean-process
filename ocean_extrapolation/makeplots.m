function [] = makeplots(plotnum)

% script to make plots for brief write up
clearvars -except plotnum; close all;
fs = 6;
cmap_bed = cptcmap('GMT_relief');

%% extended bathymetry
if plotnum==1
    
% load extended bedmachine
% (note cannot just use bedmachine because it does not cover full area)
ncfile = 'bedmachine_extended.nc';
xc = ncread(ncfile,'x');
yc = flipud(ncread(ncfile,'y'));
bed = fliplr(ncread(ncfile,'bed'));

% add bedmachine extent
bmfile = '~/Documents/BedMachineGreenland-v5.nc';
bm.x = double(ncread(bmfile,'x'));
bm.y = double(ncread(bmfile,'y'));
bm.bx = [min(bm.x),max(bm.x),max(bm.x),min(bm.x),min(bm.x)];
bm.by = [min(bm.y),min(bm.y),max(bm.y),max(bm.y),min(bm.y)];

% define ismip6 grid (in EPSG:3413)
ismip.xg = -0.72e6:1e3:0.96e6;
ismip.yg = -3.45e6:1e3:-0.57e6;
ismip.xc = 0.5*(ismip.xg(1:end-1)+ismip.xg(2:end));
ismip.yc = 0.5*(ismip.yg(1:end-1)+ismip.yg(2:end));
% ismip6 box for plotting
ismip.bx = [ismip.xg(1),ismip.xg(end),ismip.xg(end),ismip.xg(1),ismip.xg(1)];
ismip.by = [ismip.yg(1),ismip.yg(1),ismip.yg(end),ismip.yg(end),ismip.yg(1)];

t = tiledlayout(1,1,'padding','compact');
nexttile; hold on;
imagesc(xc,yc,bed'); axis xy;
colormap(cmap_bed); clim([-2000,2000]);
plot(bm.bx,bm.by,'b');
plot(ismip.bx,ismip.by,'r');
h = colorbar('fontsize',fs);
set(gca,'xtick',[],'ytick',[],'box','off');
xlim([min(xc),max(xc)]); ylim([min(yc),max(yc)]);
saveplot(7.5,10,600,'extended_bathy.png');
close all;

end

%% effective depth
if plotnum==2

% load effective depth
ncfile = 'z_eff.nc';
xc = ncread(ncfile,'x');
yc = ncread(ncfile,'y');
z_eff = double(ncread(ncfile,'z_eff'));

% read convex hull
T = readtable('convex_hull.csv');

t = tiledlayout(1,1,'padding','compact');
nexttile; hold on;
pcolor2(xc,yc,z_eff'); axis xy;
colormap(cmap_bed); clim([-2000,2000]);
plot(T.xb,T.yb,'r');
h = colorbar('fontsize',fs);
set(gca,'xtick',[],'ytick',[],'box','off','color','m');
xlim([min(xc),max(xc)]); ylim([min(yc),max(yc)]);
saveplot(7.5,10,600,'effective_depth.png');
close all;

end

%% effective positions
if plotnum==3

% load effective geometry
xy_eff_file = 'XY_eff.nc';
z_eff_file = 'z_eff.nc';
ismip.x = double(ncread(z_eff_file,'x'));
ismip.y = double(ncread(z_eff_file,'y'));
[ismip.X,ismip.Y] = ndgrid(ismip.x,ismip.y);
ismip.z_eff = double(ncread(z_eff_file,'z_eff'));
ismip.X_eff = double(ncread(xy_eff_file,'X_eff'));
ismip.Y_eff = double(ncread(xy_eff_file,'Y_eff'));
inds = find(~isnan(ismip.X_eff));
ss = 1000;

% read convex hull
T = readtable('convex_hull.csv');

t = tiledlayout(1,1,'padding','compact');
nexttile; hold on;
pcolor2(ismip.x,ismip.y,ismip.z_eff'); axis xy;
colormap(cmap_bed); clim([-2000,2000]);
plot(T.xb,T.yb,'r');
h = colorbar('fontsize',fs);
set(gca,'xtick',[],'ytick',[],'box','off','color','m');
xlim([min(ismip.x),max(ismip.x)]); ylim([min(ismip.y),max(ismip.y)]);
for i=1:ss:length(inds)
    if inpolygon(ismip.X(inds(i)),ismip.Y(inds(i)),T.xb,T.yb)
        plot([ismip.X(inds(i)),ismip.X_eff(inds(i))],[ismip.Y(inds(i)),ismip.Y_eff(inds(i))],...
            'ko-','markersize',2,'linewidth',0.25);
    end
end
saveplot(7.5,10,600,'effective_positions.png');
close all;

end

%% cmip model depths
if plotnum==4

% load effective geometry
xy_eff_file = 'XY_eff.nc';
z_eff_file = 'z_eff.nc';
ismip.x = double(ncread(z_eff_file,'x'));
ismip.y = double(ncread(z_eff_file,'y'));
[ismip.X,ismip.Y] = ndgrid(ismip.x,ismip.y);
ismip.z_eff = double(ncread(z_eff_file,'z_eff'));
ismip.X_eff = double(ncread(xy_eff_file,'X_eff'));
ismip.Y_eff = double(ncread(xy_eff_file,'Y_eff'));

% read convex hull
T = readtable('convex_hull.csv');

% load CMIP (just 1 timestep to illustrate)
cmipfile = 'thetao_Odec_MIROC6_ssp585_r1i1p1f1_gn_2020-2097.nc';
cmipfileS = 'so_Odec_MIROC6_ssp585_r1i1p1f1_gn_2020-2097.nc';
cmip.lat = ncread(cmipfile,'latitude');
cmip.lon = ncread(cmipfile,'longitude');
cmip.z0 = ncread(cmipfile,'zlev');
cmip.T0 = ncread(cmipfile,'thetao',[1,1,1,1],[Inf,Inf,Inf,1]);
cmip.S0 = ncread(cmipfileS,'so',[1,1,1,1],[Inf,Inf,Inf,1]);

% convert to EPSG:3413
[cmip.x,cmip.y] = projfwd(projcrs(3413),cmip.lat,cmip.lon);

% crop the cmip data to ISMIP ROI
xlims = [min(ismip.x),max(ismip.x)];
ylims = [min(ismip.y),max(ismip.y)];
inds = find(cmip.x>xlims(1) & cmip.x<xlims(2) & cmip.y>ylims(1) & cmip.y<ylims(2));
cmip.x = cmip.x(inds);
cmip.y = cmip.y(inds);
cmip.T0 = reshape(cmip.T0,size(cmip.T0,1)*size(cmip.T0,2),size(cmip.T0,3));
cmip.T0 = cmip.T0(inds,:);
cmip.S0 = reshape(cmip.S0,size(cmip.S0,1)*size(cmip.S0,2),size(cmip.S0,3));
cmip.S0 = cmip.S0(inds,:);

% interpolate the cmip data onto the effective depth grid
% this has to be done with care
% choose the extrap option here because we want values at the surface
% and extrapolation doesn't occur where we don't want it because
% cmip.T0 has NaNs where it doesn't have data
zvals = flipud(unique(ismip.z_eff(ismip.z_eff<=0)));
cmip.T = interp1(cmip.z0,cmip.T0',zvals,'linear','extrap')';
cmip.S = interp1(cmip.z0,cmip.S0',zvals,'linear','extrap')';

% calculate thermal forcing (without pressure/depth dependence)
l1 = -5.73e-2;
l2 = 8.32e-2;
cmip.Tfreeze_nop = l1*cmip.S+l2;
cmip.TF_nop = cmip.T - cmip.Tfreeze_nop;

% get the depth of the interpolated cmip data
for i=1:length(cmip.x)
    if sum(~isnan(cmip.TF_nop(i,:)))==0
        cmip.depth(i) = NaN;
    elseif sum(~isnan(cmip.TF_nop(i,:)))==length(zvals)
        cmip.depth(i) = zvals(end);
    else
        cmip.depth(i) = zvals(min(find(isnan(cmip.TF_nop(i,:))))-1);
    end
end
inds = find(~isnan(cmip.depth));

t = tiledlayout(1,1,'padding','compact');
nexttile; hold on;
pcolor2(ismip.x,ismip.y,ismip.z_eff'); axis xy;
scatter(cmip.x(inds),cmip.y(inds),10,cmip.depth(inds),'filled','markeredgecolor','b');
colormap(cmap_bed); clim([-2000,2000]);
plot(T.xb,T.yb,'r');
h = colorbar('fontsize',fs);
set(gca,'xtick',[],'ytick',[],'box','off','color','m');
xlim([min(ismip.x),max(ismip.x)]); ylim([min(ismip.y),max(ismip.y)]);
saveplot(7.5,10,600,'cmip_depth.png');
close all;

end

%% final result
if plotnum==5

ncfile = 'TF_example.nc';
xc = ncread(ncfile,'x');
yc = ncread(ncfile,'y');
TF = ncread(ncfile,'TF');

t = tiledlayout(1,1,'padding','compact');
nexttile; hold on;
pcolor2(xc,yc,TF'); axis xy;
h = colorbar('fontsize',fs);
set(gca,'xtick',[],'ytick',[],'box','off');
xlim([min(xc),max(xc)]); ylim([min(yc),max(yc)]);
saveplot(7.5,10,600,'TF_example.png');
close all;

end

end