clear; close all;

% load ismip bed
bedfile = 'ismip_bed.nc';
x = double(ncread(bedfile,'x'));
y = double(ncread(bedfile,'y'));
bed = double(ncread(bedfile,'bed'));
[X,Y] = ndgrid(x,y);

% load effective depth
z_eff_file = 'z_eff.nc';
z_eff = double(ncread(z_eff_file,'z_eff'));

% read convex hull
T = readtable('convex_hull.csv');

% initialise effective positions
effposition = NaN*bed;

% required masks and regions
in = double(inpolygon(X,Y,T.xb,T.yb));
belowsl = double(z_eff<=0);
in_belowsl = double(in & belowsl);
out_belowsl = double(~in & belowsl);

% split inside into connected regions
% add some artificial barriers to help delineation
bw = in_belowsl;
bw(700:900,1590)=0; % jakobshavn/peterman
% bw(479,1176:1188)=0; % jakobshavn mouth
bw(1240:1260,2249)=0; % 79N southern entry
bw(1060,820:830)=0; % tinit!
bw(1410:1435,2013)=0;
bw(1338,1710:1720)=0;
bw(1392,1720:1735)=0;
bw(1119,832)=0;
bw(860:866,420)=0;
bw(865:880,384)=0;
bw(867:870,318)=0;
bw(735,100:127)=0;
bw(709,170:175)=0;
bw(652:665,235)=0;
bw(350:358,796)=0;
bw(330,837)=0;
bw(373:378,1108)=0;
bw(420:497,1230)=0;
bw(112:160,2221)=0;
bw(371:385,1705)=0;
bw(417:430,1761)=0;
bw(407,1769)=0;
bw(363,1905:1910)=0;
bw(304,1992:2002)=0;
bw(229,2021)=0;

% do the split into regions
regions0 = bwconncomp(bw);

% version of belowsl mask that has the same modifications
bw_dist = belowsl;
inds = find(in_belowsl-bw);
bw_dist(inds) = 0;

% sort regions by size
for i=1:regions0.NumObjects
    l(i) = length(regions0.PixelIdxList{i});
end
[~,id] = sort(l,'descend');
regions.NumObjects = regions0.NumObjects;
for i=1:length(id)
    regions.PixelIdxList{i} = regions0.PixelIdxList{id(i)};
end

% get effective position for each region
% takes about 10 minutes
for k=1:regions.NumObjects
    k/regions.NumObjects
    xregion = X(regions.PixelIdxList{k});
    yregion = Y(regions.PixelIdxList{k});
    d = (xregion-mean(xregion)).^2+(yregion-mean(yregion)).^2;
    [~,id] = min(d);
    xseed = xregion(id);
    yseed = yregion(id);
    % linear index of seed point
    [~,id] = min((X(:)-xseed).^2+(Y(:)-yseed).^2);
    % connected distances from seed point
    D = bwdistgeodesic(logical(bw_dist),id,'quasi-euclidean');
    % find linear index of closest point in outside ocean
    D = D./out_belowsl;
    [~,effid] = min(D(:));
    effposition(regions.PixelIdxList{k}) = effid;
end

% add effective positions for the artificial barrier points
% by nearest neighbour interpolation
inds = find(~isnan(effposition));
f = scatteredInterpolant(X(inds),Y(inds),effposition(inds),'nearest','nearest');
inds_to_fill = find(in_belowsl-bw);
effposition(inds_to_fill) = f(X(inds_to_fill),Y(inds_to_fill));

% add effective positions for pixels outside the convex hull
effposition(find(out_belowsl)) = find(out_belowsl);

% convert from linear index to positions
X_eff = NaN*bed;
Y_eff = NaN*bed;
inds = find(~isnan(effposition));
X_eff(inds) = X(effposition(inds));
Y_eff(inds) = Y(effposition(inds));

% create netcdf
output_file = 'XY_eff.nc';
if exist(output_file,'file'), delete(output_file); end

% initialise variables in netcdf
nccreate(output_file,'x','Dimensions',{'x',length(x)},'Datatype','single');
nccreate(output_file,'y','Dimensions',{'y',length(y)},'Datatype','single');
nccreate(output_file,'X_eff','Dimensions',{'x',length(x),'y',length(y)},'Datatype','single','deflatelevel',9);
nccreate(output_file,'Y_eff','Dimensions',{'x',length(x),'y',length(y)},'Datatype','single','deflatelevel',9);

% write variables
ncwrite(output_file,'x',single(x));
ncwrite(output_file,'y',single(y));
ncwrite(output_file,'X_eff',single(X_eff));
ncwrite(output_file,'Y_eff',single(Y_eff));


% %% old code
% 
% % % do some large regions individually to save computation time
% % name = {'Petermann','Jakobshavn','Scoresby','79N','CW','NO','Humboldt','Godthab','NE','Serm','Kanger','NO2'};
% % gx = [0.0082e6,-63050,754450,476950,-299300,319450,-366800,-336800,656950,334450,544450,-33050];
% % gy = [-1.2794e6,-2280650,-1999400,-1095650,-2059400,-791900,-1058150,-2831900,-1684400,-2588150,-2336900,-791900];
% % 
% % % working plot
% % figure(); hold on;
% % pcolor2(x,y,in_belowsl');
% % axis xy equal;
% % colorbar;
% % plot(gx,gy,'ko');
% 
% for k=1:length(name)
%     % linear index of seed point
%     [~,id] = min((X(:)-gx(k)).^2+(Y(:)-gy(k)).^2);
%     % connected distances from seed point
%     D = bwdistgeodesic(logical(belowsl),id,'quasi-euclidean');
%     % find linear index of closest point in outside ocean
%     D = D./out_belowsl;
%     [~,effid] = min(D(:));
%     % find connected interior region for seed point
%     id_region = [];
%     i = 0;
%     while isempty(id_region)
%         i=i+1;
%         if any(regions.PixelIdxList{i}==id)
%             id_region = i;
%         end
%     end
%     inds = regions.PixelIdxList{id_region};
% %     % special case for peterman/jakobshavn interior regions
% %     if strcmp(name{k},'Petermann')
% %         inds(Y(inds)<-1.85e6) = [];
% %     elseif strcmp(name{k},'Jakobshavn')
% %         inds(Y(inds)>=-1.85e6) = [];
% %         inds(X(inds)<=-1.65e5) = [];
% %     end
%     % set effective position for whole region
%     effposition(inds) = effid;
% end
% 
% % the remaining effective positions have to be done individually
% remaining_inside_inds = find(belowsl & in & isnan(effposition));
% 
% tic
% for i=1:length(remaining_inside_inds)
%     i/length(remaining_inside_inds)
%     j = remaining_inside_inds(i);
%     D = bwdistgeodesic(bw,j,'quasi-euclidean');
%     D(inside_inds) = NaN;
%     if sum(~isnan(D(:)) & D(:)~=Inf)>0
%         [~,id] = min(D(:));
%         effposition(j) = id;
%     end
% end
% toc
% 
% % effective position of outside points is just the actual position
% effposition(outside_inds) = outside_inds;
% 
% %% place onto ismip6 grid
% 
% % dummy grid for now - assume ismip6 grid is the same as the grid we have
% x = bm.x;
% y = bm.y;
% X = bm.X;
% Y = bm.Y;
% 
% inds = find(~isnan(effposition(:)));
% x_eff = NaN*effposition;
% x_eff(inds) = bm.X(effposition(inds));
% y_eff = NaN*effposition;
% y_eff(inds) = bm.Y(effposition(inds));
% z_eff = effdepth;
% 
% % plot effective positions
% figure(); hold on;
% pcolor2(bm.x,bm.y,z_eff);
% for i=1:10:length(inds)
%     plot([bm.X(inds(i)),x_eff(inds(i))],[bm.Y(inds(i)),y_eff(inds(i))],'ko-');
% end
% axis xy equal; clim([-1000,0]);
% plot(T.xb,T.yb,'k');
% 
% save effective_geometry.mat x y X Y x_eff y_eff z_eff 


