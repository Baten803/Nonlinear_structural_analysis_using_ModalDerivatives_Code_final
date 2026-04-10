%% Set up data paths

baseDir = fileparts(mfilename('fullpath'));
parentDir = fullfile(baseDir, '..');

dataDir_taylor = fullfile(parentDir, ...
    'Data', 'Taylor coordinates prestress');

%% Load data

% File names
matFile_physical_allLoads = 'BridgeSection_physical_Damping_allLoads';
matFile_taylor_9_allLoads = 'BridgeSection_TaylorBasis_9_Damped_allLoads';
matFile_taylor_15_allLoads = 'BridgeSection_TaylorBasis_15_Damped_allLoads';
matFile_taylor_20_allLoads = 'BridgeSection_TaylorBasis_20_Damped_allLoads';

% Full paths
matPath_physical_allLoads = fullfile(dataDir_taylor, matFile_physical_allLoads);
matPath_taylor_9_allLoads = fullfile(dataDir_taylor, matFile_taylor_9_allLoads);
matPath_taylor_15_allLoads = fullfile(dataDir_taylor, matFile_taylor_15_allLoads);
matPath_taylor_20_allLoads = fullfile(dataDir_taylor, matFile_taylor_20_allLoads);

% Load .mat file
matData_physical_allLoads = load(matPath_physical_allLoads);
matData_taylor_9_allLoads = load(matPath_taylor_9_allLoads);
matData_taylor_15_allLoads = load(matPath_taylor_15_allLoads);
matData_taylor_20_allLoads = load(matPath_taylor_20_allLoads);

%% Full response plot 

%Setting up the physical coordinate response

T_int_phys = matData_physical_allLoads.t_CDM_phys; % Time interval
T_int_tb = matData_taylor_9_allLoads.t_CDM_tb;

T_int_phys = T_int_phys(1500000:end);
T_int_tb = T_int_tb(150000:end);

V_phys_allLoads = matData_physical_allLoads.V_phys_mid; % Nonlinear physical coordinate response 

V_phys_allLoads = V_phys_allLoads(1500000:end);


%Setting up the Taylor basis coordinates responses

%Taylor_9_allLoads
V_taylor_9_allLoads = matData_taylor_9_allLoads.V_tb_mid;
V_taylor_9_allLoads = V_taylor_9_allLoads(150000:end);
%Taylor_15_allLoads
V_taylor_15_allLoads = matData_taylor_15_allLoads.V_tb_mid;
V_taylor_15_allLoads = V_taylor_15_allLoads(150000:end);
%Taylor_20_allLoads
V_taylor_20_allLoads = matData_taylor_20_allLoads.V_tb_mid;
V_taylor_20_allLoads = V_taylor_20_allLoads(150000:end);


%% Error estimation

% Only interpolate on overlapping time interval
idx_t9 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));
idx_t15 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));
idx_t20 = (T_int_phys >= T_int_tb(1)) & (T_int_phys <= T_int_tb(end));

Vtb_9_mid_interp = nan(size(T_int_phys));
Vtb_15_mid_interp = nan(size(T_int_phys));
Vtb_20_mid_interp = nan(size(T_int_phys));

Vtb_9_mid_interp(idx_t9) = interp1(T_int_tb, V_taylor_9_allLoads, T_int_phys(idx_t9), 'linear');
Vtb_15_mid_interp(idx_t15) = interp1(T_int_tb, V_taylor_15_allLoads, T_int_phys(idx_t15), 'linear');
Vtb_20_mid_interp(idx_t20) = interp1(T_int_tb, V_taylor_20_allLoads, T_int_phys(idx_t20), 'linear');

% Error estimate for Taylor basis, only on overlap
V_phys_mid_t9 = V_phys_allLoads(idx_t9);
V_phys_mid_t15 = V_phys_allLoads(idx_t15);
V_phys_mid_t20 = V_phys_allLoads(idx_t20);

Vdiff_taylor_9 = V_phys_mid_t9 - Vtb_9_mid_interp(idx_t9);
Vdiff_taylor_15 = V_phys_mid_t15 - Vtb_15_mid_interp(idx_t15);
Vdiff_taylor_20 = V_phys_mid_t20 - Vtb_20_mid_interp(idx_t20);

e_taylor_9_MAX = max(abs(Vdiff_taylor_9)/max(abs(V_phys_mid_t9)));
e_taylor_9_AVG = (sum(abs(Vdiff_taylor_9))/length(Vdiff_taylor_9))/(sum(abs(V_phys_mid_t9))/length(V_phys_mid_t9));
e_taylor_15_MAX = max(abs(Vdiff_taylor_15)/max(abs(V_phys_mid_t15)));
e_taylor_15_AVG = (sum(abs(Vdiff_taylor_15))/length(Vdiff_taylor_15))/(sum(abs(V_phys_mid_t15))/length(V_phys_mid_t15));
e_taylor_20_MAX = max(abs(Vdiff_taylor_20)/max(abs(V_phys_mid_t20)));
e_taylor_20_AVG = (sum(abs(Vdiff_taylor_20))/length(Vdiff_taylor_20))/(sum(abs(V_phys_mid_t20))/length(V_phys_mid_t20));