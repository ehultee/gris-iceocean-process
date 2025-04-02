% script to define boundary between fjords and shelf
clear; close all;

% load ismip bed
bedfile = 'ismip_bed.nc';
x = single(ncread(bedfile,'x'));
y = single(ncread(bedfile,'y'));
b = single(ncread(bedfile,'bed')');

% make a mask for greenland and surrounding islands
t = double(b>0);
[X,Y] = meshgrid(x,y);
% remove svalbard and iceland
t(X>9e5) = 0;
% remove ellesmere island
x_elle = 1e5*[-7.4948,-6.3610,-5.0146,-4.2115,-3.7864,-3.1604,-2.7234,-1.8967,-2.3691,-4.8611,-6.3138,-7.8964];
y_elle = 1e6*[-1.1967,-1.1353,-1.0325,-0.9628,-0.9250,-0.8967,-0.8483,-0.7987,-0.6593,-0.5849,-0.5849,-0.5943];
t(inpolygon(X,Y,x_elle,y_elle)) = 0;
% remove a last bit of canada and an outlying island
t(X<-6.5e5) = 0;

% set up points for convex hull
xall = X(find(t));
yall = Y(find(t));
ss = 10;
xb = double(xall(1:ss:end));
yb = double(yall(1:ss:end));
k = boundary(xb,yb,0.5);

% plot mask and boundary
figure(); hold on;
imagesc(x,y,t); axis xy equal;
plot(xb(k),yb(k),'r','linewidth',2);

% write result
T = table(xb(k),yb(k),'VariableNames',{'xb','yb'});
writetable(T,'convex_hull.csv');

% %% old code
% 
% % make t into a mask in which 0 is water and 1 is greenland land/ice
% t(t>=3) = 0;
% t(t>1) = 1;
% % subsample dataset to make code faster (ss=1 means no subsampling)
% ss = 10;
% x = x(1:ss:end);
% y = y(1:ss:end);
% t = t(1:ss:end,1:ss:end);
% b = b(1:ss:end,1:ss:end);
% 
% % due to subsample, can get bits of water that are unconnected to the ocean
% % we want to set these to be land/ice
% not_t = 1-t;
% a = bwconncomp(not_t,8);
% for i=1:length(a.PixelIdxList),
%     l(i) = length(a.PixelIdxList{i});
% end
% [~,id] = max(l);
% inds = [];
% for i=1:length(a.PixelIdxList),
%     if l(i)~=l(id),
%         inds = [inds;a.PixelIdxList{i}];
%     end
% end
% t(inds) = 1;
% 
% % define fjord/shelf boundary and masks
% % remove islands
% % calculate connected land area (thus isolating islands)
% a = bwconncomp(t,8);
% for i=1:length(a.PixelIdxList),
%     l(i) = length(a.PixelIdxList{i});
% end
% % get mainland
% [~,id] = max(l);
% % get (x,y) coords of mainland pixels
% inds = a.PixelIdxList{id};
% [X,Y] = meshgrid(x,y);
% xb = X(inds); yb = Y(inds);
% 
% % do shrinkwrap to separate fjord and shelf
% % s=0.8 seems good for ss=10, s=0.5 for ss=5
% s = 0.3;
% k = boundary(double(xb),double(yb),s);
% xb = xb(k); yb = yb(k);
% 
% % plot to check
% figure(); hold on;
% imagesc(x,y,t); axis xy equal;
% plot(xb,yb,'r');
% 
% % write result
% writetable(table(xb,yb),'fjord_shelf_boundary.csv');


