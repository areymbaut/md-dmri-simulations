close all
clear all

%% Paths
current_dir = '/Users/alexis_reymbaut/Desktop/Research_Lund/Codes/Matrix_Gamma/Test_dtd_mv_gamma';
%odf_path = fullfile(framework_path,'tools','uvec','repulsion_angles_tri');
Watson_path = [current_dir '/Watson'];

set(0, 'defaultLegendInterpreter','latex');

%% Tested methods and metrics
% method_names = {'dtd_gamma'; 'dtd_covariance'; 'dtd_mv_gamma'};
% method_names_plot = {'Gamma'; 'Cov'; 'mv-Gamma'};
% method_color = {[1 0.65 0] [0 0.8 0] [0 0 1]}; 
% method_marker = {'o','s','h'};
method_names = {'dtd_covariance'; 'dtd_mv_gamma'};
method_names_plot = {'Covariance tensor approximation'; 'Matrix-variate Gamma approximation'};
method_color = {[1 .3 .3] [0 0 1]}; 
method_marker = {'o','h'};
param_names = {'mdiso'; 'msdanison'; 'vdiso'};
param_names_plot = {'E[$D_{\mathrm{iso}}]\;[\mu\mathrm{m}^2/\mathrm{ms}]$'; '$\tilde{\mathrm{E}}[D_{\mathrm{aniso}}^2]$'; 'V$[D_{\mathrm{iso}}]\;[\mu\mathrm{m}^4/\mathrm{ms}^2]$'};
Nparam = numel(param_names);
Nmethod = numel(method_names);

% To be able to show empty symbols for the problematic cases of dtd_covariance
method_names_overlap = {'dtd_gamma'; 'dtd_pa'; 'dtd'; 'dtd_covariance'};
method_color_overlap = {[1 0.65 0] [1 .3 .3] [0 0 1] [0 0.8 0]};
method_marker_overlap = {'o','^','*','s'};

% Only dtd_pa and dtd can access Vaniso
% method_names_Vaniso = {'dtd_pa'; 'dtd'};
% method_names_plot_Vaniso = {'MC-2D'; 'MC-4D'};
% method_color_Vaniso = {[1 .3 .3] [0 0 1]}; 
% method_marker_Vaniso = {'^','*'};
% Nmethod_Vaniso = numel(method_names_Vaniso);

%% Acquisition protocols and SNR
Nbs = 100;
acquisition_protocol_list = ["BRAIN_FWF_MERGED_mc"];
%acquisition_protocol_list = ["GE_Premier_Short", "GE_Premier_Intermediate", "GE_Premier_Long"];
SNR_list = [30, Inf];
SNR_list = sort(SNR_list);
SNR_label_list = {['=' num2str(SNR_list(1))]; '\to\infty'};
dir_protocol_list = strcat(current_dir, '/Figures_', acquisition_protocol_list);
dir_SNR_list = strcat('SNR_', string(SNR_list));
nb_SNR = length(dir_SNR_list);

%% Systems

% case_names = {'3_MeanDiso_unimodal_setvdison', '4_Viso_bimodal', '4_Viso_high_MD_bimodal', '4_Viso_unimodal_variablevdison', '5_MeanDdelta_unimodal_setvddeltan', '5_MeanDdelta_Watson_fat_structures_Ddelta', '6_VarDdelta_bimodal_ddelta', '6_VarDdelta_unimodal_variablevddeltan', '7_Dispersion_fat', '7_Dispersion_structures', '8_Crossing_two-way_real','8_Crossing_three-way_real','9_Partial_volume_real_stick_CSFf','10_Partial_volume_edema_gliomaf'};
% system_name_list = {["unimodal_setvdison_0.8","unimodal_setvdison_1.5","unimodal_setvdison_2","unimodal_setvdison_2.5"],["bimodal_1.2_1.8","bimodal_1_2","bimodal_0.7_2.3","bimodal_0.2_1.4","bimodal_0.5_2.5"],["bimodal_1.9_2.1","bimodal_1.7_2.3","bimodal_1.5_2.5","bimodal_1.4_2.6","bimodal_1.3_2.7"],["unimodal_variablevdison_0.01","unimodal_variablevdison_0.05","unimodal_variablevdison_0.1","unimodal_variablevdison_0.15"],["unimodal_setvddeltan_-0.25","unimodal_setvddeltan_0.1","unimodal_setvddeltan_0.25","unimodal_setvddeltan_0.5","unimodal_setvddeltan_0.75"],["Watson_fat_structures_Ddelta_-0.25","Watson_fat_structures_Ddelta_0.1","Watson_fat_structures_Ddelta_0.25","Watson_fat_structures_Ddelta_0.5","Watson_fat_structures_Ddelta_0.75"],["bimodal_ddelta_0.6_0.7","bimodal_ddelta_0.5_0.77","bimodal_ddelta_0.4_0.83","bimodal_ddelta_0.3_0.87","bimodal_ddelta_0.2_0.9"],["unimodal_variablevddeltan_0.01","unimodal_variablevddeltan_0.05","unimodal_variablevddeltan_0.1","unimodal_variablevddeltan_0.15"],["Watson_fat_sticks_OP_0.01","Watson_fat_sticks_OP_0.25","Watson_fat_sticks_OP_0.5","Watson_fat_sticks_OP_0.75","single_fat_stick"],["Watson_fat_structures_OP_0.01","Watson_fat_structures_OP_0.25","Watson_fat_structures_OP_0.5","Watson_fat_structures_OP_0.75","unimodal_setvddeltan_0.6816"],["two_real_sticks_crossing_15","two_real_sticks_crossing_30","two_real_sticks_crossing_45","two_real_sticks_crossing_70","two_real_sticks_crossing_90"],["three_real_sticks_crossing_15","three_real_sticks_crossing_30","three_real_sticks_crossing_45","three_real_sticks_crossing_70","three_real_sticks_crossing_90"],["composite_real_stick_and_CSFf_0","composite_real_stick_and_CSFf_0.25","composite_real_stick_and_CSFf_0.5","composite_real_stick_and_CSFf_0.75","composite_real_stick_and_CSFf_1"],["composite_edema_and_gliomaf_0","composite_edema_and_gliomaf_0.25","composite_edema_and_gliomaf_0.5","composite_edema_and_gliomaf_0.75","composite_edema_and_gliomaf_1"]};
% labels = {{'$D_{\mathrm{iso}}=0.5$','$D_{\mathrm{iso}}=1$','$D_{\mathrm{iso}}=2$','$D_{\mathrm{iso}}=2.5$'}, {'$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1,2)$','$D_{\mathrm{iso}}=(0.7,2.3)$','$D_{\mathrm{iso}}=(0.5,2.5)$'}, {'$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1,2)$','$D_{\mathrm{iso}}=(0.7,2.3)$','$D_{\mathrm{iso}}=(0.5,2.5)$'}, {'V/m2=0.01','V/m2=0.05','V/m2=0.1','V/m2=0.15'}, {'$D_\Delta=-0.25$','$D_\Delta=0.1$','$D_\Delta=0.25$','$D_\Delta=0.5$','$D_\Delta=0.75$'}, {'$D_\Delta=-0.25$','$D_\Delta=0.1$','$D_\Delta=0.25$','$D_\Delta=0.5$','$D_\Delta=0.75$'}, {'$D_\Delta=(0.6,0.7)$','$D_\Delta=(0.5,0.77)$','$D_\Delta=(0.4,0.83)$','$D_\Delta=(0.3,0.87)$','$D_\Delta=(0.2,0.9)$'}, {'V/m2=0.01','V/m2=0.05','V/m2=0.1','V/m2=0.15'}, {'$\mathrm{OP}=0.01$','$\mathrm{OP}=0.25$','$\mathrm{OP}=0.5$','$\mathrm{OP}=0.75$','$\mathrm{OP}=1$'}, {'$\mathrm{OP}=0.01$','$\mathrm{OP}=0.25$','$\mathrm{OP}=0.5$','$\mathrm{OP}=0.75$','$\mathrm{OP}=1$'}, {'$\alpha = 15^\circ$','$\alpha = 30^\circ$','$\alpha = 45^\circ$','$\alpha = 70^\circ$','$\alpha = 90^\circ$'}, {'$\alpha = 15^\circ$','$\alpha = 30^\circ$','$\alpha = 45^\circ$','$\alpha = 70^\circ$','$\alpha = 90^\circ$'}, {'$f_{\mathrm{iso}}=0$','$f_{\mathrm{iso}}=0.25$','$f_{\mathrm{iso}}=0.5$','$f_{\mathrm{iso}}=0.75$','$f_{\mathrm{iso}}=1$'}, {'$f_{\mathrm{g}}=0$','$f_{\mathrm{g}}=0.25$','$f_{\mathrm{g}}=0.5$','$f_{\mathrm{g}}=0.75$','$f_{\mathrm{g}}=1$'}};
% special_labels = {'MeanDiso', 'Viso', 'Viso', 'Viso', 'MeanDaniso', 'MeanDaniso', 'VarDaniso', 'VarDaniso', 'OP', 'OP', 'alpha2', 'alpha3', 'f_CSF', 'f_glioma'};
% x_variable_labels = {'E[$D_{\mathrm{iso}}$]'; 'V$[D_{\mathrm{iso}}]$'; 'V$[D_{\mathrm{iso}}]$'; 'V$[D_{\mathrm{iso}}]$'; '$\tilde{\mathrm{E}}[D_{\mathrm{aniso}}^2]$'; '$\tilde{\mathrm{E}}[D_{\mathrm{aniso}}^2]$'; '$\mathrm{V}_{\mathrm{A}}$'; '$\mathrm{V}_{\mathrm{A}}$'; '$\mathrm{OP}$'; '$\mathrm{OP}$'; '$\alpha$'; '$\alpha$'; '$f_{\mathrm{iso}}$'; '$f_{\mathrm{g}}$'};
case_names = {'3_MeanDiso_unimodal_setvdison', '4_Viso_bimodal', '4_Viso_high_MD_bimodal', '5_MeanDdelta_unimodal_setvddeltan', '7_Dispersion_structures','9_Partial_volume_real_stick_CSFf','10_Partial_volume_edema_gliomaf'};
system_name_list = {["unimodal_setvdison_0.8","unimodal_setvdison_1.5","unimodal_setvdison_2","unimodal_setvdison_2.5"],["bimodal_1.2_1.8","bimodal_1_2","bimodal_0.7_2.3","bimodal_0.2_1.4","bimodal_0.5_2.5"],["bimodal_1.9_2.1","bimodal_1.7_2.3","bimodal_1.5_2.5","bimodal_1.4_2.6","bimodal_1.3_2.7"],["unimodal_setvddeltan_-0.25","unimodal_setvddeltan_0.1","unimodal_setvddeltan_0.25","unimodal_setvddeltan_0.5","unimodal_setvddeltan_0.75"],["Watson_fat_structures_OP_0.01","Watson_fat_structures_OP_0.25","Watson_fat_structures_OP_0.5","Watson_fat_structures_OP_0.75","unimodal_setvddeltan_0.6816"],["composite_real_stick_and_CSFf_0","composite_real_stick_and_CSFf_0.25","composite_real_stick_and_CSFf_0.5","composite_real_stick_and_CSFf_0.75","composite_real_stick_and_CSFf_1"],["composite_edema_and_gliomaf_0","composite_edema_and_gliomaf_0.25","composite_edema_and_gliomaf_0.5","composite_edema_and_gliomaf_0.75","composite_edema_and_gliomaf_1"]};
labels = {{'$D_{\mathrm{iso}}=0.5$','$D_{\mathrm{iso}}=1$','$D_{\mathrm{iso}}=2$','$D_{\mathrm{iso}}=2.5$'}, {'$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1,2)$','$D_{\mathrm{iso}}=(0.7,2.3)$','$D_{\mathrm{iso}}=(0.5,2.5)$'}, {'$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1,2)$','$D_{\mathrm{iso}}=(0.7,2.3)$','$D_{\mathrm{iso}}=(0.5,2.5)$'}, {'$D_\Delta=-0.25$','$D_\Delta=0.1$','$D_\Delta=0.25$','$D_\Delta=0.5$','$D_\Delta=0.75$'}, {'$\mathrm{OP}=0.01$','$\mathrm{OP}=0.25$','$\mathrm{OP}=0.5$','$\mathrm{OP}=0.75$','$\mathrm{OP}=1$'}, {'$f_{\mathrm{iso}}=0$','$f_{\mathrm{iso}}=0.25$','$f_{\mathrm{iso}}=0.5$','$f_{\mathrm{iso}}=0.75$','$f_{\mathrm{iso}}=1$'}, {'$f_{\mathrm{g}}=0$','$f_{\mathrm{g}}=0.25$','$f_{\mathrm{g}}=0.5$','$f_{\mathrm{g}}=0.75$','$f_{\mathrm{g}}=1$'}};
special_labels = {'MeanDiso', 'Viso', 'Viso', 'MeanDaniso', 'OP', 'f_CSF', 'f_glioma'};
x_variable_labels = {'E[$D_{\mathrm{iso}}]\;[\mu\mathrm{m}^2/\mathrm{ms}]$'; 'V$[D_{\mathrm{iso}}]\;[\mu\mathrm{m}^4/\mathrm{ms}^2]$'; 'V$[D_{\mathrm{iso}}]\;[\mu\mathrm{m}^4/\mathrm{ms}^2]$'; '$\tilde{\mathrm{E}}[D_{\mathrm{aniso}}^2]$'; '$\mathrm{OP}$'; '$f_{\mathrm{iso}}$'; '$f_{\mathrm{g}}$'};
nb_cases = length(case_names);

%% Plotting instructions
lw = 2;
global_font_size = 13.5;
width = 0.21;
inter_h = 0.1;
height = 0.24;
x1 = (1-3*width-2*inter_h)/2;
x2 = x1 + width + inter_h;
x3 = x2 + width + inter_h;
inter_v = 0.08;

inter_error_bar = 0.06;

inter_v_top_bottom = (1 - inter_v - 2*height)/2;

y1 = inter_v_top_bottom + 0.01;
y2 = y1 + inter_v + height;
y = [y2, y1];

% For BRAIN_FWF_MERGED (shorter)
y_lim_list = [[[0.5 3],[-0.02 0.06],[-0.1 1.5]]; ... % 3_MeanDiso_unimodal_setvdison
              [[0.5 0.9000001],[-0.02 0.06],[0 0.50001]]; ... % 4_Viso_bimodal
              [[1.8 2.3000001],[-0.02 0.06],[-0.05 1.5]]; ... % 4_Viso_high_MD_bimodal
              [[0.76 0.825],[-0.05 0.8],[-0.005 0.05]]; ... % 5_MeanDdelta_unimodal_setvddeltan
              [[0.7 0.85],[0 0.65],[-0.005 0.07]]; ... % 7_Dispersion_structures
              [[0.5 3.5],[-0.01 0.55],[-0.1 5]]; ... % 9_Partial_volume_real_stick_CSFf
              [[0.5 2.75],[-0.02 0.09],[-0.1 1.75]]]; % 10_Partial_volume_edema_gliomaf

% For BRAIN_FWF_MERGED
% y_lim_list = [[[0.5 3],[-0.02 0.04],[-0.1 1.5]]; ... % 3_MeanDiso_unimodal_setvdison
%               [[0.4 1.1],[-0.05 0.15],[-0.1 2]]; ... % 4_Viso_bimodal
%               [[1.8 2.25],[-0.02 0.15],[-0.05 1.1]]; ... % 4_Viso_high_MD_bimodal
%               [[0.74 0.83],[-0.02 0.02],[-0.01 0.17]]; ... 5 4_Viso_unimodal_variablevdison
%               [[0.77 0.825],[-0.05 0.8],[-0.01 0.04]]; ... % 5_MeanDdelta_unimodal_setvddeltan
%               [[0.7 0.825],[-0.05 0.8],[-0.01 0.07]]; ... % 5_MeanDdelta_Watson_fat_structures_Ddelta
%               [[0.76 0.83],[0.35 0.55],[-0.01 0.04]]; ... % 6_VarDdelta_bimodal_ddelta
%               [[0.76 0.83],[0.4 0.55],[-0.01 0.04]]; ... % 6_VarDdelta_unimodal_variablevddeltan
%               [[0.73 0.82],[-0.01 0.60001],[-0.01 0.065]]; ... % 7_Dispersion_fat 
%               [[0.73 0.82],[-0.01 0.60001],[-0.01 0.08]]; ... % 7_Dispersion_structures
%               [[0.65 0.85],[0.1 0.50001],[-0.01 0.40001]]; ... % 8_Crossing_two-way_real
%               [[0.74 0.82],[0 0.50001],[-0.01 0.40001]]; ... % 8_Crossing_three-way_real
%               [[0.5 4],[-0.01 0.50001],[-0.1 5]]; ... % 9_Partial_volume_real_stick_CSFf
%               [[0.5 2.75],[-0.025 0.07],[-0.1 1.75]]]; % 10_Partial_volume_edema_gliomaf
          
% For in-house L6_10_16_21_S6_6_10_10_with_b0 
% y_lim_list = [[[0.5 3],[-0.075 0.25],[-0.1 1.75]]; ... % [[0.5 3],[-0.0001 0.0002],[0.005 0.06]]
%               [[0.6 1.3],[-0.1 0.5],[-0.1 3]]; ... % [[0.6 1.24],[-0.002 0.006],[0 3.1]]
%               [[1.8 2.41],[-0.05 0.15],[-0.05 2.25]]; ... % [[1.94 2.01],[-0.00051 0.0002],[-0.01 1.25]]
%               [[0.75 1.1],[-0.1 0.35],[-0.01 1.5]]; ...
%               [[0.7 1.1],[-0.1 1],[-0.1 1]]; ... % [[0.75 1.2],[-0.1 1],[-0.1 1.5]]
%               [[0.72 1],[-0.1 1],[-0.1 1]]; ... % [[0.77 0.81],[-0.1 1],[-0.04 0.01]]
%               [[0.75 1.1],[0.2 0.65],[-0.1 1]]; ...
%               [[0.75 1.05],[0.2 0.70001],[-0.1 1]]; ...
%               [[0.73 1],[0.25 0.75],[-0.1 1]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.73 1],[0.2 0.75],[-0.1 1]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.72 1],[0.15 0.6],[-0.025 1]]; ...
%               [[0.72 1.1],[0.2 0.61],[-0.1 1]]; ...
%               [[0.5 2.75],[-0.05 0.55],[-0.05 2.5]]; ...
%               [[0.5 2.75],[-0.075 0.3],[-0.05 3]]];

% For in-house_no_reg_for_cov (and Rician noise)
% y_lim_list = [[[0.5 3.5],[-0.05 0.2000001],[-0.1 3.5]]; ... % [[0.5 3],[-0.0001 0.0002],[0.005 0.06]]
%               [[0.65 1.200001],[-0.05 0.40000001],[-0.1 3]]; ... % [[0.6 1.24],[-0.002 0.006],[0 3.1]]
%               [[1.85 2.45],[-0.05 0.40000001],[-0.1 3]]; ... % [[1.94 2.01],[-0.00051 0.0002],[-0.01 1.25]]
%               [[0.75 1],[-0.05 0.3],[-0.01 1.25]]; ...
%               [[0.75 0.95],[-0.05 0.85],[-0.05 0.65]]; ... % [[0.75 1.2],[-0.1 1],[-0.1 1.5]]
%               [[0.75 0.95],[-0.05 0.85],[-0.05 0.65]]; ... % [[0.77 0.81],[-0.1 1],[-0.04 0.01]]
%               [[0.75 1.1],[0.2 0.65],[-0.1 1]]; ...
%               [[0.75 1.05],[0.2 0.70001],[-0.1 1]]; ...
%               [[0.74 1],[0.2 0.7],[-0.1 0.75]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.74 0.96],[0.3 0.75],[-0.05 0.7]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.74 1],[0.2 0.600001],[-0.1 1]]; ...
%               [[0.74 1],[0.2 0.600001],[-0.1 1]]; ...
%               [[0.5 4],[-0.05 0.5],[-0.1 5.5]]; ...
%               [[0.5 2.75],[-0.05 0.3],[-0.1 2]]];

% For in-house_no_reg_for_cov (and Rician noise and do_weight for Gamma)
% y_lim_list = [[[0.5 3.1],[-0.05 0.2000001],[-0.1 1.75]]; ... % [[0.5 3],[-0.0001 0.0002],[0.005 0.06]]
%               [[0.65 1.200001],[-0.05 0.40000001],[-0.1 2.5]]; ... % [[0.6 1.24],[-0.002 0.006],[0 3.1]]
%               [[1.85 2.400001],[-0.05 0.40000001],[-0.1 2.5]]; ... % [[1.94 2.01],[-0.00051 0.0002],[-0.01 1.25]]
%               [[0.75 1],[-0.05 0.3],[-0.01 1.25]]; ...
%               [[0.75 1],[-0.05 0.85],[-0.05 0.7]]; ... % [[0.75 1.2],[-0.1 1],[-0.1 1.5]]
%               [[0.75 1],[-0.05 0.85],[-0.05 0.7]]; ... % [[0.77 0.81],[-0.1 1],[-0.04 0.01]]
%               [[0.75 1.1],[0.2 0.65],[-0.1 1]]; ...
%               [[0.75 1.05],[0.2 0.70001],[-0.1 1]]; ...
%               [[0.74 1],[0.2 0.7],[-0.1 0.75]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.74 1.00001],[0.3 0.700001],[-0.05 0.800001]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.74 1],[0.2 0.600001],[-0.1 1]]; ...
%               [[0.74 1],[0.2 0.600001],[-0.1 1]]; ...
%               [[0.5 3.75],[-0.05 0.500001],[-0.1 3]]; ...
%               [[0.5 2.75],[-0.05 0.3],[-0.1 2]]];

% For in-house_no_reg_for_cov
% y_lim_list = [[[0.5 3],[-0.05 0.2000001],[-0.1 1.5]]; ... % [[0.5 3],[-0.0001 0.0002],[0.005 0.06]]
%               [[0.7 1.200001],[-0.05 0.45],[-0.1 3]]; ... % [[0.6 1.24],[-0.002 0.006],[0 3.1]]
%               [[1.85 2.3],[-0.05 0.45],[-0.1 3]]; ... % [[1.94 2.01],[-0.00051 0.0002],[-0.01 1.25]]
%               [[0.75 1],[-0.05 0.3],[-0.01 1.25]]; ...
%               [[0.7 1],[-0.05 0.8],[-0.1 0.75]]; ... % [[0.75 1.2],[-0.1 1],[-0.1 1.5]]
%               [[0.7 1],[-0.05 0.8],[-0.1 0.75]]; ... % [[0.77 0.81],[-0.1 1],[-0.04 0.01]]
%               [[0.75 1.1],[0.2 0.65],[-0.1 1]]; ...
%               [[0.75 1.05],[0.2 0.70001],[-0.1 1]]; ...
%               [[0.74 1],[0.2 0.7],[-0.1 0.75]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.74 0.95],[0.3 0.7],[-0.1 0.6]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.74 1],[0.2 0.600001],[-0.1 1]]; ...
%               [[0.74 1],[0.2 0.600001],[-0.1 1]]; ...
%               [[0.5 3.5],[-0.05 0.5],[-0.1 4]]; ...
%               [[0.5 2.75],[-0.05 0.3],[-0.1 2]]];

% For in-house_with_reg_for_cov
% y_lim_list = [[[0.5 3],[-0.025 0.2],[-0.1 0.8]]; ... % [[0.5 3],[-0.0001 0.0002],[0.005 0.06]]
%               [[0.7 1.200001],[-0.05 0.4],[-0.1 3]]; ... % [[0.6 1.24],[-0.002 0.006],[0 3.1]]
%               [[1.6 2.41],[-0.05 0.4],[-0.1 3]]; ... % [[1.94 2.01],[-0.00051 0.0002],[-0.01 1.25]]
%               [[0.75 1],[-0.05 0.3],[-0.01 1.25]]; ...
%               [[0.75 0.975],[-0.05 0.8],[-0.1 0.7]]; ... % [[0.75 1.2],[-0.1 1],[-0.1 1.5]]
%               [[0.75 0.975],[-0.05 0.8],[-0.1 0.7]]; ... % [[0.77 0.81],[-0.1 1],[-0.04 0.01]]
%               [[0.75 1.1],[0.2 0.65],[-0.1 1]]; ...
%               [[0.75 1.05],[0.2 0.70001],[-0.1 1]]; ...
%               [[0.73 1],[0.2 0.7],[-0.1 0.75]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.73 1],[0.3 0.7],[-0.1 0.75]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.74 1],[0.2 0.600001],[-0.1 1]]; ...
%               [[0.74 1],[0.2 0.600001],[-0.1 1]]; ...
%               [[0.5 2.75],[-0.05 0.55],[-0.1 1.75]]; ...
%               [[0.5 2.75],[-0.05 0.3],[-0.1 2]]];
          
% For Intermediate GE
% y_lim_list = [[[0.5 3],[-0.075 0.25],[-0.1 1.75]]; ... % [[0.5 3],[-0.0001 0.0002],[0.005 0.06]]
%               [[0.6 1.24],[-0.1 0.5],[0 3.1]]; ... % [[0.6 1.24],[-0.002 0.006],[0 3.1]]
%               [[1.75 2.5],[-0.05 0.175],[-0.01 1.25]]; ... % [[1.94 2.01],[-0.00051 0.0002],[-0.01 1.25]]
%               [[0.7 1.1],[-0.05 0.35],[-0.01 0.75]]; ...
%               [[0.75 1.2],[-0.1 1],[-0.1 1.5]]; ... % [[0.75 1.2],[-0.1 1],[-0.1 1.5]]
%               [[0.72 1.2],[-0.1 1],[-0.1 1.5]]; ... % [[0.77 0.81],[-0.1 1],[-0.04 0.01]]
%               [[0.75 1.2],[0.1 0.7],[-0.1 1.75]]; ...
%               [[0.75 1.25],[0.1 0.7],[-0.1 1.75]]; ...
%               [[0.73 1.22],[0.15 0.75],[-0.1 1.75]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.73 1.12],[0.15 0.75],[-0.1 1.5]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.7 0.975],[0.2 0.8],[-0.025 0.3]]; ...
%               [[0.74 0.975],[0.2 0.65],[-0.025 0.42]]; ...
%               [[0.5 2.75],[-0.05 0.55],[-0.05 2.7]]; ...
%               [[0.5 2.75],[-0.075 0.3],[-0.05 2.2]]];

% For LP
% y_lim_list = [[[0.5 3.5],[-0.075 0.4],[-0.1 1.75]]; ... % [[0.5 3],[-0.0001 0.0002],[0.005 0.06]]
%               [[0.6 1.4],[-0.1 0.6],[0 4]]; ... % [[0.6 1.24],[-0.002 0.006],[0 3.1]]
%               [[1.75 2.5],[-0.05 0.2],[-0.01 2]]; ... % [[1.94 2.01],[-0.00051 0.0002],[-0.01 1.25]]
%               [[0.7 1],[-0.1 0.4],[-0.01 0.75]]; ...
%               [[0.75 1],[-0.1 1],[-0.1 0.5]]; ... % [[0.75 1.2],[-0.1 1],[-0.1 1.5]]
%               [[0.72 1],[-0.1 1],[-0.1 0.5]]; ... % [[0.77 0.81],[-0.1 1],[-0.04 0.01]]
%               [[0.75 1],[0.1 0.7],[-0.1 0.5]]; ...
%               [[0.75 1],[0.1 0.7],[-0.1 0.5]]; ...
%               [[0.73 1],[0.15 0.75],[-0.1 0.5]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.73 1.05],[0.15 0.85],[-0.1 0.5]]; ... % [[0.78 0.805],[0.15 0.75],[-0.02 0.005]]
%               [[0.7 1],[0.1 0.65],[-0.025 0.3]]; ...
%               [[0.7 1.1],[0.1 0.65],[-0.025 0.5]]; ...
%               [[0.5 2.75],[-0.15 0.7],[-0.05 2.7]]; ...
%               [[0.5 2.75],[-0.075 0.4],[-0.05 3]]];

error_cap_method = [10 8 6 4];
          
%% Main
for it_acq_prot = 1:length(dir_protocol_list)
    dir_protocol = dir_protocol_list(it_acq_prot);
    for it_case = 1:nb_cases
        case_name = case_names{it_case}
        nb_systems = length(labels{it_case});
        if ~strcmp(special_labels{it_case},'')
            dir = strcat(dir_protocol, '/', dir_SNR_list(1));
            matfile_names = strcat(dir, '/', case_name, '_SNR', num2str(SNR_list(1)), '_system', string(1:nb_systems), '.mat');
            labels_system = cell([1, nb_systems]);
            if strcmp(special_labels{it_case},'MeanDiso') || strcmp(special_labels{it_case},'MeanDaniso') || strcmp(special_labels{it_case},'Viso') || strcmp(special_labels{it_case},'VarDaniso')
                for it_system = 1:nb_systems    
                    structure_metrics = load(matfile_names(it_system));
                    if strcmp(special_labels{it_case},'MeanDiso')
                        new_label_system = structure_metrics.mdiso_true;
                        labels_system{it_system} = ['$' num2str(round(new_label_system*1e9,2)) '$'];
                    end
                    
                    if strcmp(special_labels{it_case},'MeanDaniso')
                        new_label_system = structure_metrics.msdanison_true;
                        labels_system{it_system} = ['$' num2str(round(new_label_system,2)) '$'];
                    end
                    
                    if strcmp(special_labels{it_case},'Viso')
                        new_label_system = structure_metrics.vdiso_true*1e18;
                        labels_system{it_system} = ['$' num2str(round(new_label_system,2)) '$'];
                    end
                    
                    if strcmp(special_labels{it_case},'VarDaniso')
                        new_label_system = structure_metrics.vsdaniso_true*1e36;
                        labels_system{it_system} = ['$' num2str(round(new_label_system,2)) '$'];
                    end
                end
            elseif strcmp(special_labels{it_case},'OP')
                OP_list = ["$0.01$", "$0.25$", "$0.5$", "$0.75$", "$1$"];
                for it_system = 1:nb_systems
                    labels_system{it_system} = char(OP_list(it_system));
                end
            elseif strcmp(special_labels{it_case},'alpha2') || strcmp(special_labels{it_case},'alpha3')
                angle_list = ["$15$", "$30$", "$45$", "$75$", "$90$"];
                for it_system = 1:nb_systems
                    labels_system{it_system} = char(angle_list(it_system));
                end
            elseif strcmp(special_labels{it_case},'f_CSF') || strcmp(special_labels{it_case},'f_glioma')
                fraction_list = ["$0$", "$0.25$", "$0.5$", "$0.75$", "$1$"];
                for it_system = 1:nb_systems
                    labels_system{it_system} = char(fraction_list(it_system));
                end
            end
        else
            labels_system = labels{it_case};
        end
        
        zeroarray = zeros(nb_systems,nb_SNR);
        for nparam = 1:Nparam
            for nmethod = 1:Nmethod
                eval(['median_' param_names{nparam} '_' method_names{nmethod} ' = zeroarray;'])
                eval(['Q1_' param_names{nparam} '_' method_names{nmethod} ' = zeroarray;'])
                eval(['Q3_' param_names{nparam} '_' method_names{nmethod} ' = zeroarray;'])
            end
        end
        
        clear x_variable
        x_variable = zeros([1, nb_systems]);
        
        for nparam = 1:Nparam
            eval(['clear true_' param_names{nparam}])
        end
            
        for it_system = 1:nb_systems
            x_variable(it_system) = str2double(extractBetween(labels_system{it_system},'$','$'));
            for it_SNR = 1:nb_SNR
                dir = strcat(dir_protocol, '/', dir_SNR_list(it_SNR));
                matfile = strcat(dir, '/', case_name, '_SNR', num2str(SNR_list(it_SNR)), '_system', num2str(it_system), '.mat');
                structure_metrics = load(matfile);
                for nparam = 1:Nparam
                    if strcmp(param_names{nparam},'vsdaniso') || strcmp(param_names{nparam},'vsdanison')
%                         eval(['true = structure_metrics.', param_names{nparam}, '_true;'])
%                         eval(['true_' param_names{nparam} '(it_system, it_SNR) = true;'])
%                         for nmethod = 1:Nmethod_Vaniso
%                             eval(['quartiles = quantile(structure_metrics.', param_names{nparam}, '_', method_names_Vaniso{nmethod}, ',3);'])
%                             median = quartiles(2);
%                             Q1 = quartiles(1);
%                             Q3 = quartiles(3);
%                             eval(['median_' param_names{nparam} '_' method_names_Vaniso{nmethod} '(it_system, it_SNR) = median;'])
%                             eval(['Q1_' param_names{nparam} '_' method_names_Vaniso{nmethod} '(it_system, it_SNR) = Q1;'])
%                             eval(['Q3_' param_names{nparam} '_' method_names_Vaniso{nmethod} '(it_system, it_SNR) = Q3;'])
%                         end
                    else
                        eval(['true = structure_metrics.', param_names{nparam}, '_true;'])
                        eval(['true_' param_names{nparam} '(it_system, it_SNR) = true;'])
                        for nmethod = 1:Nmethod
                            eval(['quartiles = quantile(structure_metrics.', param_names{nparam}, '_', method_names{nmethod}, ',3);'])
                            median = quartiles(2);
                            Q1 = quartiles(1);
                            Q3 = quartiles(3);
                            eval(['median_' param_names{nparam} '_' method_names{nmethod} '(it_system, it_SNR) = median;'])
                            eval(['Q1_' param_names{nparam} '_' method_names{nmethod} '(it_system, it_SNR) = Q1;'])
                            eval(['Q3_' param_names{nparam} '_' method_names{nmethod} '(it_system, it_SNR) = Q3;'])
                        end
                    end
                end
            end
        end
        
        if length(unique(x_variable)) < length(x_variable)
            x_variable(1) = -x_variable(1);
        end
        
        numeric_labels_system = str2double(extractBetween(labels_system,'$','$'));
        if length(unique(numeric_labels_system)) < length(numeric_labels_system)
            numeric_labels_system(1) = -numeric_labels_system(1);
            txt = extractBetween(labels_system{1},'$','$');
            labels_system{1} = ['$' num2str(round(-str2double(txt{1}),2)) '$'];
        end
        
        f = figure('Position', [0,0,1100,1500], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        set(gcf, 'Renderer', 'painters')
        
        for it_SNR = 1:2
       
            axh_mdiso = axes('position',[x1 y(it_SNR) width height]);
            axh_msdanison = axes('position',[x2 y(it_SNR) width height]);
            axh_vdiso = axes('position',[x3 y(it_SNR) width height]);
            
            clear L
            
            for nparam = 1:Nparam
                eval(['axh = axh_' param_names{nparam} ';'])
                hold(axh,'on')
                
                axh.XAxis.TickLabelInterpreter = 'latex';
                axh.YAxis.TickLabelInterpreter = 'latex';
                
                if strcmp(param_names{nparam},'mdiso')
                    correct_unit = 1e9;
                elseif strcmp(param_names{nparam},'vdiso')
                    correct_unit = 1e18;
                elseif strcmp(param_names{nparam},'vsdaniso')
                    correct_unit = 1e36;
                else
                    correct_unit = 1;
                end
                
                eval(['true_to_plot = true_' param_names{nparam} '(:, it_SNR);'])
                plot(axh, x_variable, correct_unit*true_to_plot, '--', 'Linewidth', 3*lw/4, 'Color', [0, 0, 0]) %[0.75, 0, 0.75]
                
                %             if strcmp(param_names{nparam},'vsdaniso') || strcmp(param_names{nparam},'vsdanison')
                %                 for nmethod = 1:Nmethod_Vaniso
                %                     if strcmp(method_names_Vaniso{nmethod}, 'dtd')
                %                         for it_system = 1:nb_systems
                %                             eval(strcat('plot(axh, it_system+delta, correct_unit*', stat, '_', param_names{nparam}, '_', method_names_Vaniso{nmethod}, '(it_system, :), ', strcat('''-', string(method_marker_Vaniso{nmethod}), ''''), ', ', '''MarkerSize''', ', ', string(6.55), ', ', '''MarkerFaceColor''', ', ', string('method_color_Vaniso{nmethod}'), ', ', '''Linewidth''', ', ', num2str(lw/2), ', ', '''Color''', ', ', string('method_color_Vaniso{nmethod}'), ')'))
                %                         end
                %                     else
                %                         for it_system = 1:nb_systems
                %                             eval(strcat('plot(axh, it_system+delta, correct_unit*', stat, '_', param_names{nparam}, '_', method_names_Vaniso{nmethod}, '(it_system, :), ', strcat('''-', string(method_marker_Vaniso{nmethod}), ''''), ', ', '''MarkerFaceColor''', ', ', string('method_color_Vaniso{nmethod}'), ', ', '''Linewidth''', ', ', num2str(lw/2), ', ', '''Color''', ', ', string('method_color_Vaniso{nmethod}'), ')'))
                %                         end
                %                     end
                %                 end
                %             else
                %                 for nmethod = 1:Nmethod
                %                     if strcmp(method_names_overlap{nmethod}, 'dtd')
                %                         for it_system = 1:nb_systems
                %                             eval(strcat('plot(axh, it_system+delta, correct_unit*', stat, '_', param_names{nparam}, '_', method_names_overlap{nmethod}, '(it_system, :), ', strcat('''-', string(method_marker_overlap{nmethod}), ''''), ', ', '''MarkerSize''', ', ', string(6.55), ', ', '''MarkerFaceColor''', ', ', string('method_color_overlap{nmethod}'), ', ', '''Linewidth''', ', ', num2str(lw/2), ', ', '''Color''', ', ', string('method_color_overlap{nmethod}'), ')'))
                %                         end
                %                     else
                %                         for it_system = 1:nb_systems
                %                             eval(strcat('plot(axh, it_system+delta, correct_unit*', stat, '_', param_names{nparam}, '_', method_names_overlap{nmethod}, '(it_system, :), ', strcat('''-', string(method_marker_overlap{nmethod}), ''''), ', ', '''MarkerFaceColor''', ', ', string('method_color_overlap{nmethod}'), ', ', '''Linewidth''', ', ', num2str(lw/2), ', ', '''Color''', ', ', string('method_color_overlap{nmethod}'), ')'))
                %                         end
                %                     end
                %                 end
                %                 L(nmethod) = plot(axh, NaN,1, ['-' method_marker{nmethod}],'color',method_color{nmethod},'MarkerFaceColor',method_color{nmethod},'MarkerEdgeColor',method_color{nmethod}); %// dummy plot for legend
                %             end
                
                for nmethod = 1:Nmethod
                    for it_system = 1:nb_systems
                        eval(strcat('medians = median_', param_names{nparam}, '_', method_names{nmethod}, '(:, it_SNR);'))
                        eval(strcat('Q1s = Q1_', param_names{nparam}, '_', method_names{nmethod}, '(:, it_SNR);'))
                        eval(strcat('Q3s = Q3_', param_names{nparam}, '_', method_names{nmethod}, '(:, it_SNR);'))
                        lower = medians - Q1s;
                        upper = Q3s - medians;
                        if it_SNR == 1
%                             eval(strcat('p1= plot(axh, x_variable, correct_unit*Q1s, ', '''Color''', ', ', string('method_color{nmethod}'), ', ', '''Linewidth''', ', ', num2str(lw/2), ');'))
%                             eval(strcat('p2= plot(axh, x_variable, correct_unit*Q3s, ', '''Color''', ', ', string('method_color{nmethod}'), ', ', '''Linewidth''', ', ', num2str(lw/2), ');'))
                            eval(strcat('p3= plot(axh, x_variable, correct_unit*medians, ', '''Color''', ', ', string('method_color{nmethod}'), ', ', '''Linewidth''', ', ', num2str(lw/4), ');'))
%                             p1.Color(4) = 0.1;
%                             p2.Color(4) = 0.1;
                            p3.Color(4) = 0;
                            
                            eval(strcat('errorbar(axh, x_variable, correct_unit*medians, zeros(size(lower)), zeros(size(upper)), ', strcat('''', string(method_marker{nmethod}), ''''), ', ', '''MarkerFaceColor''', ', ', string('method_color{nmethod}'), ', ', '''Capsize''', ', ', num2str(0.1), ', ', '''Linewidth''', ', ', num2str(3*lw/4), ', ', '''Color''', ', ', string('method_color{nmethod}'), ')'))
                            eval(strcat('errorbar(axh, x_variable+inter_error_bar*(max(x_variable)-min(x_variable)), correct_unit*medians, correct_unit*lower, correct_unit*upper, ', strcat('''', ' ', ''''), ', ', '''MarkerFaceColor''', ', ', string('method_color{nmethod}'), ', ', '''Capsize''', ', ', num2str(error_cap_method(nmethod)), ', ', '''Linestyle''', ', ', string('''none'''), ', ', '''Linewidth''', ', ', num2str(3*lw/4), ', ', '''Color''', ', ', string('method_color{nmethod}'), ')'))

%                             eval(strcat('errorbar(axh, x_variable, correct_unit*medians, correct_unit*lower, correct_unit*upper, ', strcat('''', string(method_marker{nmethod}), ''''), ', ', '''MarkerFaceColor''', ', ', string('method_color{nmethod}'), ', ', '''Linewidth''', ', ', num2str(3*lw/4), ', ', '''Color''', ', ', string('method_color{nmethod}'), ')'))
                        else
                            eval(strcat('errorbar(axh, x_variable, correct_unit*medians, zeros(size(lower)), zeros(size(upper)), ', strcat('''', string(method_marker{nmethod}), ''''), ', ', '''MarkerFaceColor''', ', ', string('method_color{nmethod}'), ', ', '''Capsize''', ', ', num2str(0), ', ', '''Linewidth''', ', ', num2str(3*lw/4), ', ', '''Color''', ', ', string('method_color{nmethod}'), ')'))
                        end
                    end
                    L(nmethod) = plot(axh, NaN,1, ['-' method_marker{nmethod}],'color',method_color{nmethod},'MarkerFaceColor',method_color{nmethod},'MarkerEdgeColor',method_color{nmethod}); %// dummy plot for legend
                end
                
                box(axh,'off')
                grid(axh,'on')
                
                if it_SNR == 2
                    set(axh,'XTick',numeric_labels_system)
                    set(axh,'XTickLabel',labels_system)
                    xlabel(axh,x_variable_labels{it_case}, 'Fontsize', global_font_size, 'Interpreter', 'latex')
                else
                    set(axh,'XTick',numeric_labels_system)
                    set(axh,'XTickLabel',{' '})
                end
                
                set(axh, 'Fontsize', global_font_size)
                xlim(axh, [x_variable(1) x_variable(end)+inter_error_bar*(max(x_variable)-min(x_variable))])
                
%                 if it_SNR == 1
%                     ylim(axh, [y_lim_list(it_case, 2*nparam-1) y_lim_list(it_case, 2*nparam)])
%                 end
                
                ylim(axh, [y_lim_list(it_case, 2*nparam-1) y_lim_list(it_case, 2*nparam)])
                
                if it_SNR == 1
                    ylab = ylabel(axh,param_names_plot{nparam}, 'Fontsize', global_font_size, 'Interpreter', 'latex');
                    set(ylab, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
                end
                
                if nparam == 1
                   tle = title(axh, ['$\mathrm{SNR} ' SNR_label_list{it_SNR} '$'], 'Fontsize', global_font_size, 'Interpreter', 'latex');
                   set(tle, 'Units', 'Normalized', 'Position', [0.5, 1.05, 0]);
                end
                
                set(axh, 'LineWidth', lw)
                set(axh, 'TickDir','out');
                if nparam == 1 && it_SNR == 2
                    [h_legend, hObj] = legend(axh, L, method_names_plot);
                    pos = get(h_legend,'Position');
                    set(h_legend,'Position', [x1, y(2) - 4*height/6, pos(3:4)]);
                    hl = findobj(hObj,'type','line');  
                    set(hl,'LineWidth',1.5);
                end
                
            end
        end
        saveas(gcf,strcat(dir_protocol, '/metrics_', case_name, '.pdf'))
        clf(f)
    end
end