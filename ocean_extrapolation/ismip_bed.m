% script to get bathymetry on ismip7 grid
clear; close all;

% define ismip7 grid (in EPSG:3413)
ismip.res = 1e3;
% grid cell edges
ismip.xg = -0.72e6:1e3:0.96e6;
ismip.yg = -3.45e6:1e3:-0.57e6;
% grid cell centres
ismip.xc = 0.5*(ismip.xg(1:end-1)+ismip.xg(2:end));
ismip.yc = 0.5*(ismip.yg(1:end-1)+ismip.yg(2:end));
% ismip7 box for plotting
ismip.bx = [ismip.xg(1),ismip.xg(end),ismip.xg(end),ismip.xg(1),ismip.xg(1)];
ismip.by = [ismip.yg(1),ismip.yg(1),ismip.yg(end),ismip.yg(end),ismip.yg(1)];
% initialise bed
ismip.bed = NaN(length(ismip.xc),length(ismip.yc));

% load extended bedmachine
% (note cannot just use bedmachine because it does not cover full area)
ncfile = 'bedmachine_extended.nc';
bm.xc = ncread(ncfile,'x');
bm.yc = flipud(ncread(ncfile,'y'));
bm.bed = fliplr(ncread(ncfile,'bed'));
% note that we have done some flipping and transposing so that the bed
% is bed(x,y) and y is an increasing vector
bm.res = 150;
% create a vector of grid cell edges for extended bedmachine
bm.xg = [bm.xc-0.5*bm.res;bm.xc(end)+0.5*bm.res];
bm.yg = [bm.yc-0.5*bm.res;bm.yc(end)+0.5*bm.res];

% plot to illustrate
figure(); hold on;
imagesc(bm.xc,bm.yc,bm.bed'); axis xy equal;
plot(ismip.bx,ismip.by,'k');

% get bathymetry on ismip7 grid by block-averaging
% takes around 7 minutes
ss = 1;
tic
for i=1:ss:length(ismip.xc)
    i/length(ismip.xc)
    for j=1:ss:length(ismip.yc)

        % find high-res grid points for this ismip7 grid cell
        i_inds = find(bm.xc>ismip.xg(i)-0.5*bm.res & bm.xc<ismip.xg(i+1)+0.5*bm.res);
        j_inds = find(bm.yc>ismip.yg(j)-0.5*bm.res & bm.yc<ismip.yg(j+1)+0.5*bm.res);

        % calculate weight matrix
        wx = ones(length(i_inds),1);
        wx(1) = (bm.xg(i_inds(2))-ismip.xg(i))/bm.res;
        wx(end) = (ismip.xg(i+1)-bm.xg(i_inds(end)))/bm.res;
        wy = ones(1,length(j_inds));
        wy(1) = (bm.yg(j_inds(2))-ismip.yg(j))/bm.res;
        wy(end) = (ismip.yg(j+1)-bm.yg(j_inds(end)))/bm.res;
        w = (wx.*wy)*bm.res^2/ismip.res^2;

        % calculate weighted bed value
        ismip.bed(i,j) = sum(bm.bed(i_inds,j_inds).*w,'all');

    end
end
toc

% plot to check
figure(); hold on;
pcolor2(ismip.xc,ismip.yc,ismip.bed'); axis xy equal;
h = colorbar;

% create netcdf
output_file = 'ismip_bed.nc';
if exist(output_file, 'file'), delete(output_file); end

% initialise variables in netcdf
nccreate(output_file,'x','Dimensions',{'x',length(ismip.xc)},'Datatype','single');
nccreate(output_file,'y','Dimensions',{'y',length(ismip.yc)},'Datatype','single');
nccreate(output_file,'bed','Dimensions',{'x',length(ismip.xc),'y',length(ismip.yc)},'Datatype','single','deflatelevel',9);

% write variables
ncwrite(output_file,'x',single(ismip.xc));
ncwrite(output_file,'y',single(ismip.yc));
ncwrite(output_file,'bed',single(ismip.bed));
