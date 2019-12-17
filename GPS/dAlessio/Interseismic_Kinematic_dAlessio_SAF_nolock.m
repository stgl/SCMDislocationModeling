%%% Initialize program %%%
disp(' ')
disp('%%% BEM Kinematic SCM, dAlessio SAF, BAVU GPS Stations %%%')
disp(' ')
t_tot=0;    
 c=clock;                %comment before sending to CEES
 PID=feature('getpid');  %comment before sending to CEES
disp(['Date/time in Stanford= ',num2str(c(1)),' ',num2str(c(2)),' ',...
    num2str(c(3)),' ',num2str(c(4)),' ',num2str(c(5)),' ',num2str(c(6))])
disp(['PID= ',num2str(PID)])


addpath(genpath('/data/cees/hilley/SCM_Aron/Gmsh'))
% addpath(genpath('/home/faron/bin/BEM/tribemx'))
addpath(genpath('/data/cees/aron/bin/matlab'))

% addpath(genpath('/Users/felipearon/Documents/POSTDOC/Codes/BEM/Jacks_tribemx/tribem2018'))
% addpath(genpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/Gmsh'))
% addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/CEES_runs')
% addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/Gmsh/dAlessio')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Loading files%%%%

tic

% %%%%%% comment before sending to CEES
 run_tag='a';
% % dAlessio strike-slips in km/yr. Negative is dextral. Same order as faultnames vector.
 v=3;
 G_out=0;
 gps_file='obs_data_gps.txt';
 grid_file='obs_data_grid.txt';
 save_output=1;
% %%%%%%

xyz_gps=load(gps_file);
xyz_grid=load(grid_file);
xyz=[xyz_gps;xyz_grid];
len_gps=num2str(length(xyz_gps));
len_grid=num2str(length(xyz_grid));
%     disp('grids loaded= chi_obs && grid_obs')
disp(['Grids loaded= ',gps_file,' & ',grid_file])
gridsname=[gps_file(1:end-4),'_',grid_file(1:end-4)];

disp(['"v" parameter output structure= ',num2str(v)])

obs_data=struct('x',xyz(:,1),'y',xyz(:,2),'z',xyz(:,3),'v',v);
filename = 'SAF_nolock_kinematic_model_out';

x_min = 521.95; 
x_max = 633.82;
y_min = 4072.72; 
y_max = 4206.88;

dx = 0.5;

[X,Y] = meshgrid([x_min:dx:x_max],[y_min:dx:y_max]);
sn = size(X);
nx = sn(1);
ny = sn(2);

xyz = [reshape(X,prod(size(X)),1) reshape(Y,prod(size(Y)),1) zeros(prod(size(Y)),1)];

obs_data=struct('x',xyz(:,1),'y',xyz(:,2),'z',xyz(:,3),'v',v);

% Faults long term interseismic slip
SAF_filename= 'SAF_SF_Pen_SCM_SJB_nolock.msh';

faultnames=[SAF_filename];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Patching faults input structure%%%

faults=ReadPatches(faultnames);

ends = cumsum(faults.nEl); % Give ending indices of the n faults
begs = [1; ends(1:end-1)+1]; % Give beginning indices of the n faults
centroids=PatchCentroid(faults.c,faults.v); % xyz of faults centroids, km ENU

disp(['Nr of elem= ',num2str(sum(faults.nEl))])
% disp(['Nr of free parameters= ',num2str(2*(sum(faults.nEl(1:end-nr_box_patches))))])
disp('Nr of free parameters= 0')
disp(['Nr obs points= ',num2str(length(xyz)),' (gps= ',len_gps,', grid= ',len_grid,')'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assigning kinematic BCs dAlessio model %%%% 

d = zeros(sum(faults.nEl), 3); % Initialize the boundary condition array
bc = d;
d = [17.1 0 0.2;
     17.1 0 0.2;
     16.8 0 3.3;
     16.8 0 3.3;
     16.4 0 4.9;
     16.4 0 4.9;
     16.4 0 4.9;
     16.4 0 4.9;
     16.9 0 2.5;
     16.9 0 2.5;
     16.9 0 2.5;
     16.9 0 2.5;
     16.9 0 2.5;
     16.9 0 2.5;
     16.9 0 2.5;
     16.9 0 2.5;
     16.9 0 2.5;
     16.9 0 2.5;
     20.2 0 -1.9;
     20.2 0 -1.9;] * -1E-6;

d(:,3) = 0.0;

%%%Running kinematic interseismic model%%%%
% tic
disp('Run kinematic interseismic dAlessio model...')

if G_out
    disp('Saving G ic-matrix to workspace= true')
%    [slip_kine, trac_kine, obs_kine, G] = tribemx(faults, d, bc, obs_data);   
else
    disp('Saving G ic-matrix to workspace= false')
    [slip_kine, trac_kine, obs_kine] = tribemx(faults, d, bc, obs_data);
end

%%%Saving files%%%%
% tic

if save_output
%    filename = [run_tag,'_output_SCM_kinematic_dAlessio_',num2str(c(1)),'_',num2str(c(2)),'_',...
    num2str(c(3));
    filename = strcat(filename,['_PID',num2str(PID),'.mat']);    
%    save(filename,'slip_kine','trac_kine','obs_kine')
    % save(['G_cl',cl,'_chi_grid.mat'],'G','-v7.3')
    disp('Saving outputs= true')
else
    disp('Saving outputs= false')
end

t=toc;
t_tot=t_tot+t;
% disp(['saving shear Greens Fcns and G matrix takes ',num2str(t),' secs'])
disp(['Running interseismic kinematic model and saving outputs takes ',num2str(t),' secs'])

in_use = monitor_memory_whos;

u = obs_kine.u;
uz = reshape(obs_kine.u,3,length(u)/3)';
UZ = reshape(uz(:,3),nx,ny);


disp(['Memory in workspace after running shear model= ',num2str(in_use),' MB'])

disp(['Program takes a total of ',num2str(t_tot),' secs'])
disp('good bye!')
