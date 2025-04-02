clear; close all;

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

% for each ismip point find the closest deep enough point
ismip.x_ind = NaN*ismip.X;
ismip.z_ind = NaN*ismip.X;
% loop over effective depths
for k=1:length(zvals)
    k
    % get all points with this effective depth
    i_inds = find(ismip.z_eff==zvals(k));
    % save the effective depth index
    ismip.z_ind(i_inds) = k;
    % get all cmip points that have data at this depth
    c_inds = find(cmip.depth<=zvals(k));
    % find closest eligible cmip point for each ismip point
    for i=1:length(i_inds)
        % distance from cmip points to this ismip point
        d = (cmip.x(c_inds)-ismip.X_eff(i_inds(i))).^2+(cmip.y(c_inds)-ismip.Y_eff(i_inds(i))).^2;
        % closest one
        [~,id] = min(d);
        % save index of closest cmip point
        ismip.x_ind(i_inds(i)) = c_inds(id);
    end
end

% now it follows that TF on the ismip points are given by
% ismip.TF_nop = cmip.TF_nop(ismip.x_ind,ismip.z_ind)
% probably can speed this up by vectorising
% but it doesn't take long
ismip.TF_nop = NaN*ismip.X;
for i=1:length(ismip.x)
    for j=1:length(ismip.y)
        if ~isnan(ismip.x_ind(i,j))
            ismip.TF_nop(i,j) = cmip.TF_nop(ismip.x_ind(i,j),ismip.z_ind(i,j));
        end
    end
end

% finally correct for the freezing point dependence on pressure
l3 = 7.61e-4;
ismip.TF = ismip.TF_nop - l3*ismip.z_eff;
% check if any are less than 0 and fix
inds = find(ismip.TF<0);
ismip.TF(inds) = 0;
if length(inds)>0
    disp(['Warning: ',num2str(length(inds)),' pixels corrected for TF<0']);
end

% lastly, all the unconnected below sea level points must be assigned
% TF=0. These points have z_eff=NaN
ismip.TF(isnan(ismip.z_eff))=0;

figure();
pcolor2(ismip.x,ismip.y,ismip.TF');
axis xy equal; h=colorbar;
ylabel(h,'thermal forcing (C)');

% create netcdf
output_file = 'TF_example.nc';
if exist(output_file, 'file'), delete(output_file); end

% initialise variables in netcdf
nccreate(output_file,'x','Dimensions',{'x',length(ismip.x)},'Datatype','single');
nccreate(output_file,'y','Dimensions',{'y',length(ismip.y)},'Datatype','single');
nccreate(output_file,'TF','Dimensions',{'x',length(ismip.x),'y',length(ismip.y)},'Datatype','single','deflatelevel',9);

% write variables
ncwrite(output_file,'x',single(ismip.x));
ncwrite(output_file,'y',single(ismip.y));
ncwrite(output_file,'TF',single(ismip.TF));

% % plot to visualise
% figure(); hold on;
% pcolor2(ismip.x,ismip.y,ismip.z_eff');
% plot_inds = find(cmip.depth<=0);
% plot(cmip.x(plot_inds),cmip.y(plot_inds),'ko');
% xplot = [ismip.X(i_inds(i)),ismip.X_eff(i_inds(i)),cmip.x(ismip.ind(i_inds(i)))];
% yplot = [ismip.Y(i_inds(i)),ismip.Y_eff(i_inds(i)),cmip.y(ismip.ind(i_inds(i)))];
% plot(xplot,yplot,'ro-');
% plot(T.xb,T.yb,'m');
% axis equal;
% colorbar;

