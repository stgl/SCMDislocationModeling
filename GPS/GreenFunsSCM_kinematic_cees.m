%%% Initialize program %%%
disp(' ')
disp('%%% BEM Kinematic SCM BAVU GPS Stations Greens %%%')
disp(' ')
t_tot=0;    
% c=clock;                %comment before sending to CEES
% PID=feature('getpid');  %comment before sending to CEES
disp(['Date/time in Stanford= ',num2str(c(1)),' ',num2str(c(2)),' ',...
    num2str(c(3)),' ',num2str(c(4)),' ',num2str(c(5)),' ',num2str(c(6))])
disp(['PID= ',num2str(PID)])


addpath(genpath('/data/cees/aron/SCM/Gmsh'))
% addpath(genpath('/home/faron/bin/BEM/tribemx'))
addpath(genpath('/data/cees/aron/bin/matlab'))

% addpath(genpath('/Users/felipearon/Documents/POSTDOC/Codes/BEM/Jacks_tribemx/tribem2018'))
% addpath(genpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/Gmsh'))
% addpath('/Users/felipearon/Documents/POSTDOC/BEM Santa Cruz Mnts/tribemx_tries/CEES_runs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Loading files%%%%

tic

% %%%%%% comment before sending to CEES
% cl='2.5';
% locking_depth=-10;  %km, ENU system
% creeping_N=4075;    %km, North end coordinate creping section SA, ENU system
% v=3;
% G_out=0;
% gps_file='obs_data_gps.txt';
% grid_file='obs_data_grid.txt';
% greens_shear_longterm_file='GreenFncs_cl2.5_shear.mat';
% greens_push_longterm_file='GreenFncs_cl2.5_push.mat';
% save_output=1;
% %%%%%%

disp(['Fault elements characteristic length= ',cl,' km'])
disp(['Faults locking depth= ',num2str(-locking_depth),' km'])
if creeping_N
    disp(['SA creeping Northern end coordinate= ',num2str(creeping_N),' km (UTM)'])
end

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


%%% Faults
BE_filename= ['BE_cl',cl,'.msh'];
SHAMTV_filename= ['SHAMTV_cl',cl,'.msh'];
SA_filename= ['SA_cl',cl,'.msh'];

% Decollement
NAdecoll_filename= ['NA_decoll_cl',cl,'.msh'];
PAdecoll_filename= ['PA_decoll_cl',cl,'.msh'];

% Box patches (tectonic BCs)
BoxL_filename= ['box_L_cl',cl,'.msh'];
BoxR_filename= ['box_R_cl',cl,'.msh'];
BoxLL_filename= ['box_LL_cl',cl,'.msh'];
BoxUR_filename= ['box_UR_cl',cl,'.msh'];
BoxUL_filename= ['box_UL_cl',cl,'.msh'];
BoxLR_filename= ['box_LR_cl',cl,'.msh'];

nr_faults=3;
nr_decoll=2;
nr_box_patches=6;

faultnames=[BE_filename,' ',SHAMTV_filename,' ',SA_filename, ' ',...
    NAdecoll_filename, ' ',PAdecoll_filename,' ',BoxL_filename,' ',...
    BoxR_filename,' ',BoxLL_filename,' ',BoxUR_filename,' ',...
    BoxUL_filename,' ',BoxLR_filename];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Patching faults input structure%%%
faults=ReadPatches(faultnames);

ends = cumsum(faults.nEl); % Give ending indices of the 3 faults
begs = [1; ends(1:end-1)+1]; % Give beginning indices of the 3 faults

centroids=PatchCentroid(faults.c,faults.v); % xyz of faults centroids, km ENU

bc = zeros(sum(faults.nEl), 3); % Initialize the boundary condition type array.
%Interseismic kinematic model does not have free parameters so bc=[0 0 0]='bbb'

disp(['Nr of elem= ',num2str(sum(faults.nEl))])
% disp(['Nr of free parameters= ',num2str(2*(sum(faults.nEl(1:end-nr_box_patches))))])
disp('Nr of free parameters= 0')
disp(['Nr obs points= ',num2str(length(xyz)),' (gps= ',len_gps,', grid= ',len_grid,')'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Assigning kinematic shear tectonic BCs%%%% 

load(greens_shear_longterm_file);

d = slip_shear; % Initialize the boundary condition array
% Lock patches on faults above the locking depth

d(centroids(begs(1):ends(3),3)>=locking_depth & centroids(begs(1):ends(3),2)>=creeping_N,:)=...
    d(centroids(begs(1):ends(3),3)>=locking_depth & centroids(begs(1):ends(3),2)>=creeping_N,:)*0;

%%%Running shear model%%%%
% tic
disp('Run kinematic shear model...')

if G_out
    disp('Saving G ic-matrix to workspace= true')
    [slip_shear_kine, trac_shear_kine, obs_shear_kine, G] = tribemx(faults, d, bc, obs_data);   
else
    disp('Saving G ic-matrix to workspace= false')
    [slip_shear_kine, trac_shear_kine, obs_shear_kine] = tribemx(faults, d, bc, obs_data);
end

%%%Saving files%%%%
% tic

if save_output
    filename = [run_tag,'_output_SCM_Greens_shear_kinematic_',num2str(c(1)),'_',num2str(c(2)),'_',...
    num2str(c(3))];
    filename = strcat(filename,['_PID',num2str(PID),'.mat']);    
    save(filename,'obs_shear_kine')
    % save(['G_cl',cl,'_chi_grid.mat'],'G','-v7.3')
    disp('Saving Green Fcns= true')
else
    disp('Saving Green Fcns= false')
end

t=toc;
t_tot=t_tot+t;
% disp(['saving shear Greens Fcns and G matrix takes ',num2str(t),' secs'])
disp(['Running shear kinematic model and saving Greens Funs takes ',num2str(t),' secs'])

in_use = monitor_memory_whos;
disp(['Memory in workspace after running shear model= ',num2str(in_use),' MB'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Assigning kinematic push tectonic BCs%%%% 
tic
load(greens_push_longterm_file);

d = slip_push; % Initialize the boundary condition array
% Lock patches on faults above the locking depth

d(centroids(begs(1):ends(3),3)>=locking_depth,:)=...
    d(centroids(begs(1):ends(3),3)>=locking_depth,:)*0;

%%%Running push model%%%%
% tic
disp('Run kinematic push model...')

if G_out
    disp('Saving G ic-matrix to workspace= true')
    [slip_push_kine, trac_push_kine, obs_push_kine, G] = tribemx(faults, d, bc, obs_data);   
else
    disp('Saving G ic-matrix to workspace= false')
    [slip_push_kine, trac_push_kine, obs_push_kine] = tribemx(faults, d, bc, obs_data);
end

%%%Saving files%%%%
% tic

if save_output
    filename = [run_tag,'_output_SCM_Greens_push_kinematic_',num2str(c(1)),'_',num2str(c(2)),'_',...
    num2str(c(3))];
    filename = strcat(filename,['_PID',num2str(PID),'.mat']);    
    save(filename,'obs_push_kine')
    % save(['G_cl',cl,'_chi_grid.mat'],'G','-v7.3')
    disp('Saving Green Fcns= true')
else
    disp('Saving Green Fcns= false')
end
t=toc;
t_tot=t_tot+t;
% disp(['saving shear Greens Fcns and G matrix takes ',num2str(t),' secs'])
disp(['Running push kinematic model and saving Greens Funs takes ',num2str(t),' secs'])

in_use = monitor_memory_whos;
disp(['Memory in workspace after running push model= ',num2str(in_use),' MB'])

disp(['Program takes a total of ',num2str(t_tot),' secs'])
disp('good bye!')
