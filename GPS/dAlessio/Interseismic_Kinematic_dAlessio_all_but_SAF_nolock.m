addpath(genpath('../../Gmsh'))
addpath(genpath('/data/cees/aron/bin/matlab'))

slip_creep_vec=[2.7, 6.6, 4.0, 6.2, 12.9, 6.9, 5.4, 5.4].*-1e-6;

%gps_file='obs_data_gps.txt';

%xyz_gps=load(gps_file);

v=3;

x_min = 500;
x_max = 680;
y_min = 4040;
y_max = 4220;

dx = 0.5;
[X,Y] = meshgrid([x_min:dx:x_max],[y_min:dx:y_max]);
sn = size(X);
nx = sn(1);
ny = sn(2);
xyz = [reshape(X,prod(size(X)),1) reshape(Y,prod(size(Y)),1) zeros(prod(size(Y)),1)];

obs_data=struct('x',xyz(:,1),'y',xyz(:,2),'z',xyz(:,3),'v',v);

faultnames=['SG_slip_nolock.msh', ' ', ...
  'RC_H_slip_nolock.msh', ' ', 'WN_slip_nolock.msh', ' ', 'C_north_slip_nolock.msh', ' ', ...
  'C_cen_south_slip_nolock.msh', ' ', 'GV_Con_slip_nolock.msh', ' ', 'Gr_slip_nolock.msh', ' ', ...
  'VM_slip_nolock.msh'];

faults=ReadPatches(faultnames);

ends = cumsum(faults.nEl); % Give ending indices of the n faults
begs = [1; ends(1:end-1)+1]; % Give beginning indices of the n faults
centroids=PatchCentroid(faults.c,faults.v); % xyz of faults centroids, km ENU
bc = zeros(sum(faults.nEl), 3);

d = zeros(sum(faults.nEl), 3);
for i=1:length(faults.nEl)
    d(begs(i):ends(i), :) = repmat([slip_creep_vec(i) 0 0], faults.nEl(i), 1); % Assign [strike-slip 0 0] as the magnitude of boundary conditions to faults
end

d(:,3) = 0.0;

[slip_kine, trac_kine, obs_kine] = tribemx(faults, d, bc, obs_data);

u = reshape(obs_kine.u,3,length(obs_kine.u)/3)';
UX = reshape(u(:,1),nx,ny);
UY = reshape(u(:,2),nx,ny);
UZ = reshape(u(:,3),nx,ny);

meshview(faults.c, faults.v, d(:,1));
hold on
maxval = max(abs(UZ(:)));
imagesc(X(1,:),Y(:,1), UZ, [-maxval maxval]);
colormap jet
colorbar

figure
quiver(X,Y,UX,UY);
axis equal;

figure
imagesc(X(1,:),Y(:,1),UZ,[-maxval maxval]);
colormap jet
colorbar
set(gca,'ydir','normal');axis equal;axis image
