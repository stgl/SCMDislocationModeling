%%% Initialize program %%%
disp(' ')
disp(' ')
t_tot=0;    
 c=clock;                %comment before sending to CEES
 PID=feature('getpid');  %comment before sending to CEES
disp(['Date/time in Stanford= ',num2str(c(1)),' ',num2str(c(2)),' ',...
    num2str(c(3)),' ',num2str(c(4)),' ',num2str(c(5)),' ',num2str(c(6))])
disp(['PID= ',num2str(PID)])


addpath(genpath('/data/cees/hilley/SCM_Aron/Gmsh'))
addpath(genpath('/data/cees/aron/bin/matlab'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Loading files%%%%

tic

% %%%%%% comment before sending to CEES
 run_tag='a';
% % dAlessio strike-slips in km/yr. Negative is dextral. Same order as faultnames vector.
 v=3;
 G_out=0;
 save_output=1;
% %%%%%%

filename = 'Straight_out';

x_min = -20.00; 
x_max = 20.00;
y_min = -20.00; 
y_max = 20.00;

dx = 0.1;
[X,Y] = meshgrid([x_min:dx:x_max],[y_min:dx:y_max]);
sn = size(X);
nx = sn(1);
ny = sn(2);

xyz = [reshape(X,prod(size(X)),1) reshape(Y,prod(size(Y)),1) zeros(prod(size(Y)),1)];

obs_data=struct('x',xyz(:,1),'y',xyz(:,2),'z',xyz(:,3),'v',v);

% Faults long term interseismic slip
SAF_filename= 'StraightDislocation.msh';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assigning kinematic BCs dAlessio model %%%% 

d = zeros(sum(faults.nEl), 3); % Initialize the boundary condition array
bc = d;
d = [1.0 0 0.0;
     1.0 0 0.0;] * -1;

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
us = reshape(obs_kine.u,3,length(u)/3)';
UX = reshape(us(:,1),nx,ny);
UY = reshape(us(:,2),nx,ny);
UZ = reshape(us(:,3),nx,ny);


disp(['Memory in workspace after running shear model= ',num2str(in_use),' MB'])

disp(['Program takes a total of ',num2str(t_tot),' secs'])
disp('good bye!')
