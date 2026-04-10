%% Set up data paths

clear all
close all

baseDir = fileparts(mfilename('fullpath'));
parentDir = fullfile(baseDir, '..');

dataDir_taylor = fullfile(parentDir, ...
    'Data', 'Taylor coordinates with V0');
% C:\Users\matth\OneDrive\Skrivebord\Kandidat projekt\Multi dof system scripts\Beam systems\

%% Load data

% File names
matFile_physical_withV0 = 'BridgeSection_physical_withV0_damped_allLoads';
matFile_taylor_3_withV0 = 'BridgeSection_TaylorBasis_3_withV0_Damping_allLoads';
matFile_taylor_6_withV0 = 'BridgeSection_TaylorBasis_6_withV0_Damping_allLoads';
matFile_taylor_9_withV0 = 'BridgeSection_TaylorBasis_9_withV0_Damping_allLoads';
matFile_taylor_3_withV0_atV0 = 'BridgeSection_TaylorBasis_3_withV0_atV0_Damping_allLoads';
matFile_taylor_6_withV0_atV0 = 'BridgeSection_TaylorBasis_6_withV0_atV0_Damping_allLoads';
matFile_taylor_9_withV0_atV0 = 'BridgeSection_TaylorBasis_9_withV0_atV0_Damping_allLoads';

% Full paths
matPath_physical_withV0 = fullfile(dataDir_taylor, matFile_physical_withV0);
matPath_taylor_3_withV0 = fullfile(dataDir_taylor, matFile_taylor_3_withV0);
matPath_taylor_6_withV0 = fullfile(dataDir_taylor, matFile_taylor_6_withV0);
matPath_taylor_9_withV0 = fullfile(dataDir_taylor, matFile_taylor_9_withV0);
matPath_taylor_3_withV0_atV0 = fullfile(dataDir_taylor, matFile_taylor_3_withV0_atV0);
matPath_taylor_6_withV0_atV0 = fullfile(dataDir_taylor, matFile_taylor_6_withV0_atV0);
matPath_taylor_9_withV0_atV0 = fullfile(dataDir_taylor, matFile_taylor_9_withV0_atV0);

% Load .mat file
matData_physical_withV0 = load(matPath_physical_withV0);
matData_taylor_3_withV0 = load(matPath_taylor_3_withV0);
matData_taylor_6_withV0 = load(matPath_taylor_6_withV0);
matData_taylor_9_withV0 = load(matPath_taylor_9_withV0);
matData_taylor_3_withV0_atV0 = load(matPath_taylor_3_withV0_atV0);
matData_taylor_6_withV0_atV0 = load(matPath_taylor_6_withV0_atV0);
matData_taylor_9_withV0_atV0 = load(matPath_taylor_9_withV0_atV0);

%% Full response plot

%Setting up the physical coordinate response

T_int_phys = matData_physical_withV0.t_CDM_phys; % Time interval
T_int_tb = matData_taylor_3_withV0.t_CDM_tb;

V_phys_withV0 = matData_physical_withV0.V_phys_mid; % Nonlinear physical coordinate response 


%Setting up the Taylor basis coordinates responses

%Taylor_3_withV0
V_taylor_3_withV0 = matData_taylor_3_withV0.VsT_mid;
%Taylor_6_withV0
V_taylor_6_withV0 = matData_taylor_6_withV0.VsT_mid;     
%Taylor_9_withV0
V_taylor_9_withV0 = matData_taylor_9_withV0.VsT_mid;
%Taylor_3_withhhV0_atV0
V_taylor_3_withV0_atV0 = matData_taylor_3_withV0_atV0.VsT_mid;
%Taylor_6_withV0_atV0
V_taylor_6_withV0_atV0 = matData_taylor_6_withV0_atV0.VsT_mid;
%Taylor_9_withV0_atV0
V_taylor_9_withV0_atV0 = matData_taylor_9_withV0_atV0.VsT_mid;

%% Error estimation for V = 0

% Only interpolate on overlapping time interval
idx_t3_withV0 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));
idx_t6_withV0 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));
idx_t9_withV0 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));

Vtb_3_mid_interp_withV0 = nan(size(T_int_phys));
Vtb_6_mid_interp_withV0 = nan(size(T_int_phys));
Vtb_9_mid_interp_withV0 = nan(size(T_int_phys));

Vtb_3_mid_interp_withV0(idx_t3_withV0) = interp1(T_int_tb, V_taylor_3_withV0, T_int_phys(idx_t3_withV0), 'linear');
Vtb_6_mid_interp_withV0(idx_t6_withV0) = interp1(T_int_tb, V_taylor_6_withV0, T_int_phys(idx_t6_withV0), 'linear');
Vtb_9_mid_interp_withV0(idx_t9_withV0) = interp1(T_int_tb, V_taylor_9_withV0, T_int_phys(idx_t9_withV0), 'linear');

% Error estimate for Taylor basis, only on overlap
V_phys_mid_t3_withV0 = V_phys_withV0(idx_t3_withV0);
V_phys_mid_t6_withV0 = V_phys_withV0(idx_t6_withV0);
V_phys_mid_t9_withV0 = V_phys_withV0(idx_t9_withV0);

Vdiff_taylor_3_withV0 = V_phys_mid_t3_withV0 - Vtb_3_mid_interp_withV0(idx_t3_withV0);
Vdiff_taylor_6_withV0 = V_phys_mid_t6_withV0 - Vtb_6_mid_interp_withV0(idx_t6_withV0);
Vdiff_taylor_9_withV0 = V_phys_mid_t9_withV0 - Vtb_9_mid_interp_withV0(idx_t9_withV0);

e_taylor_3_withV0_MAX = max(abs(Vdiff_taylor_3_withV0)/max(abs(V_phys_mid_t3_withV0)));
e_taylor_3_withV0_AVG = (sum(abs(Vdiff_taylor_3_withV0))/length(Vdiff_taylor_3_withV0))/(sum(abs(V_phys_mid_t3_withV0))/length(V_phys_mid_t3_withV0));
e_taylor_6_withV0_MAX = max(abs(Vdiff_taylor_6_withV0)/max(abs(V_phys_mid_t6_withV0)));
e_taylor_6_withV0_AVG = (sum(abs(Vdiff_taylor_6_withV0))/length(Vdiff_taylor_6_withV0))/(sum(abs(V_phys_mid_t6_withV0))/length(V_phys_mid_t6_withV0));
e_taylor_9_withV0_MAX = max(abs(Vdiff_taylor_9_withV0)/max(abs(V_phys_mid_t9_withV0)));
e_taylor_9_withV0_AVG = (sum(abs(Vdiff_taylor_9_withV0))/length(Vdiff_taylor_9_withV0))/(sum(abs(V_phys_mid_t9_withV0))/length(V_phys_mid_t9_withV0));

%% Error estimation for V = V0

% Only interpolate on overlapping time interval
idx_t3_withV0_atV0 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));
idx_t6_withV0_atV0 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));
idx_t9_withV0_atV0 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));

Vtb_3_mid_interp_withV0_atV0 = nan(size(T_int_phys));
Vtb_6_mid_interp_withV0_atV0 = nan(size(T_int_phys));
Vtb_9_mid_interp_withV0_atV0 = nan(size(T_int_phys));

Vtb_3_mid_interp_withV0_atV0(idx_t3_withV0_atV0) = interp1(T_int_tb, V_taylor_3_withV0_atV0, T_int_phys(idx_t3_withV0_atV0), 'linear');
Vtb_6_mid_interp_withV0_atV0(idx_t6_withV0_atV0) = interp1(T_int_tb, V_taylor_6_withV0_atV0, T_int_phys(idx_t6_withV0_atV0), 'linear');
Vtb_9_mid_interp_withV0_atV0(idx_t9_withV0_atV0) = interp1(T_int_tb, V_taylor_9_withV0_atV0, T_int_phys(idx_t9_withV0_atV0), 'linear');

% Error estimate for Taylor basis, only on overlap
V_phys_mid_t3_withV0_atV0 = V_phys_withV0(idx_t3_withV0_atV0);
V_phys_mid_t6_withV0_atV0 = V_phys_withV0(idx_t6_withV0_atV0);
V_phys_mid_t9_withV0_atV0 = V_phys_withV0(idx_t9_withV0_atV0);

Vdiff_taylor_3_withV0_atV0 = V_phys_mid_t3_withV0_atV0 - Vtb_3_mid_interp_withV0_atV0(idx_t3_withV0_atV0);
Vdiff_taylor_6_withV0_atV0 = V_phys_mid_t6_withV0_atV0 - Vtb_6_mid_interp_withV0_atV0(idx_t6_withV0_atV0);
Vdiff_taylor_9_withV0_atV0 = V_phys_mid_t9_withV0_atV0 - Vtb_9_mid_interp_withV0_atV0(idx_t9_withV0_atV0);

e_taylor_3_withV0_atV0_MAX = max(abs(Vdiff_taylor_3_withV0_atV0)/max(abs(V_phys_mid_t3_withV0_atV0)));
e_taylor_3_withV0_atV0_AVG = (sum(abs(Vdiff_taylor_3_withV0_atV0))/length(Vdiff_taylor_3_withV0_atV0))/(sum(abs(V_phys_mid_t3_withV0_atV0))/length(V_phys_mid_t3_withV0_atV0));
e_taylor_6_withV0_atV0_MAX = max(abs(Vdiff_taylor_6_withV0_atV0)/max(abs(V_phys_mid_t6_withV0_atV0)));
e_taylor_6_withV0_atV0_AVG = (sum(abs(Vdiff_taylor_6_withV0_atV0))/length(Vdiff_taylor_6_withV0_atV0))/(sum(abs(V_phys_mid_t6_withV0_atV0))/length(V_phys_mid_t6_withV0_atV0));
e_taylor_9_withV0_atV0_MAX = max(abs(Vdiff_taylor_9_withV0_atV0)/max(abs(V_phys_mid_t9_withV0_atV0)));
e_taylor_9_withV0_atV0_AVG = (sum(abs(Vdiff_taylor_9_withV0_atV0))/length(Vdiff_taylor_9_withV0_atV0))/(sum(abs(V_phys_mid_t9_withV0_atV0))/length(V_phys_mid_t9_withV0_atV0));