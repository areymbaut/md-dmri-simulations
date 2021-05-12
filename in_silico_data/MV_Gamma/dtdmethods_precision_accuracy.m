close all
clear all

%% Paths
current_dir = '/Users/alexis_reymbaut/Desktop/Research_Lund/Codes/Matrix_Gamma/Test_dtd_mv_gamma';
framework_path = '/Users/alexis_reymbaut/Desktop/Research_Lund/Codes/md-dmri/md-dmri-master';
pa_path = fullfile(framework_path,'tools','uvec','repulsion_angles');
Watson_path = [current_dir '/Watson'];
%odf_path = fullfile(framework_path,'tools','uvec','repulsion_angles_tri');

%% Tested methods and metrics
method_names = {'dtd_gamma'; 'dtd_covariance'; 'dtd_mv_gamma'};
method_names_plot = {'Gamma'; 'Cov'; 'mv-Gamma'};
method_color = {[1 0.65 0] [0 0.8 0] [0 0 1]}; 
param_names = {'mdiso'; 'msdanison'; 'vdiso'};
param_names_plot = {'E$_{\mathrm{I}}$'; '$\tilde{\mathrm{E}}_{\mathrm{A}}$'; '$\mathrm{V}_{\mathrm{I}}$'};
Nparam = numel(param_names);
Nmethod = numel(method_names);

mkdir(current_dir, 'Distributions_systems');

Nbs = 100; % Number of inversion realizations
SNR_list = [Inf];

%% Loop on acquisition protocols

acquisition_protocol_list = ["BRAIN_FWF_MERGED_mc"];
%acquisition_protocol_list = ["GE_Premier_Short", "GE_Premier_Intermediate", "GE_Premier_Long"];
% acquisition_protocol_list = ["Sherbrooke", "Pseudo_Rand1", "Pseudo_Rand2", "Pseudo_Rand3", "Pseudo_Rand4", "Pseudo_Rand5", "Pseudo_Rand6", "Pseudo_Rand7", "Pseudo_Rand8"]; 
% 'all_15', '10_10_15_20', 'Sherbrooke', 'Pseudo_Rand'

for it_acq_prot = 1:length(acquisition_protocol_list)
    
    acquisition_protocol = acquisition_protocol_list(it_acq_prot);
%     xps = load([current_dir '/in_house_xps.mat']);
%     
%     xps = load([current_dir '/xps_b0_6_10_21_46_extremes_mix_b0_lin_0_10_16_21_sph_2_2_2_0_pla_0_0_16_21.mat']);
%     xps = xps.xps;
    
    xps = mdm_xps_load('BRAIN_FWF_MERGED_mc_xps.mat');
    
    mkdir(current_dir, strcat('Figures_', acquisition_protocol));
    figure_dir_parent = strcat(current_dir, '/Figures_', acquisition_protocol);
    
    acquisition_protocol
    
    %% Loop on SNRs
    for SNR = SNR_list
        
        SNR
        
        mkdir(figure_dir_parent, strcat('SNR_', num2str(SNR)));
        figure_dir = strcat(figure_dir_parent, '/SNR_', num2str(SNR));
        
        timestart = tic;
        
        %% Mean[Diso] system
%         system_name_list = ["unimodal_setvdison_0.8","unimodal_setvdison_1.5","unimodal_setvdison_2","unimodal_setvdison_2.5"];
%         labels_system = {'$D_{\mathrm{iso}}=0.8$','$D_{\mathrm{iso}}=1.5$','$D_{\mathrm{iso}}=2$','$D_{\mathrm{iso}}=2.5$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, '', labels_system, figure_dir, ['/3_MeanDiso_boxplot_unimodal_setvdison_SNR' num2str(SNR) '.pdf'], method_color, 'MeanDiso')
% % %         
% % %         %% Var[Diso] systems
%         system_name_list = ["bimodal_0.7_0.9","bimodal_0.5_1.1","bimodal_0.3_1.3","bimodal_0.2_1.4","bimodal_0.1_1.5"];
%         labels_system = {'$D_{\mathrm{iso}}=(1.2,1.8)$','$D_{\mathrm{iso}}=(1,2)$','$D_{\mathrm{iso}}=(0.7,2.3)$','$D_{\mathrm{iso}}=(0.7,2.3)$','$D_{\mathrm{iso}}=(0.5,2.5)$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, '', labels_system, figure_dir, ['/4_Viso_boxplot_bimodal_SNR' num2str(SNR) '.pdf'], method_color, 'Viso')
% %         
%         system_name_list = ["bimodal_1.9_2.1","bimodal_1.7_2.3","bimodal_1.5_2.5","bimodal_1.4_2.6","bimodal_1.3_2.7"];
%         labels_system = {'$D_{\mathrm{iso}}=(1.9,2.1)$','$D_{\mathrm{iso}}=(1.7,2.3)$','$D_{\mathrm{iso}}=(1.5,2.5)$','$D_{\mathrm{iso}}=(1.4,2.6)$','$D_{\mathrm{iso}}=(1.3,2.7)$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, '', labels_system, figure_dir, ['/4_Viso_high_MD_boxplot_bimodal_SNR' num2str(SNR) '.pdf'], method_color, 'Viso')
%         
%         system_name_list = ["unimodal_variablevdison_0.01","unimodal_variablevdison_0.05","unimodal_variablevdison_0.1","unimodal_variablevdison_0.15"];
%         labels_system = {'V/m2=0.01','V/m2=0.05','V/m2=0.1','V/m2=0.15'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, '', labels_system, figure_dir, ['/4_Viso_boxplot_unimodal_variablevdison_SNR' num2str(SNR) '.pdf'], method_color, 'Viso')
% %         
% %         %% Mean[Ddelta] systems
%         system_name_list = ["unimodal_setvddeltan_-0.25","unimodal_setvddeltan_0.1","unimodal_setvddeltan_0.25","unimodal_setvddeltan_0.5","unimodal_setvddeltan_0.75"];
%         labels_system = {'$D_\Delta=-0.25$','$D_\Delta=0.1$','$D_\Delta=0.25$','$D_\Delta=0.5$','$D_\Delta=0.75$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, '', labels_system, figure_dir, ['/5_MeanDdelta_boxplot_unimodal_setvddeltan_SNR' num2str(SNR) '.pdf'], method_color, 'MeanDaniso')
%          
%         system_name_list = ["Watson_fat_structures_Ddelta_-0.25","Watson_fat_structures_Ddelta_0.1","Watson_fat_structures_Ddelta_0.25","Watson_fat_structures_Ddelta_0.5","Watson_fat_structures_Ddelta_0.75"];
%         labels_system = {'$D_\Delta=-0.25$','$D_\Delta=0.1$','$D_\Delta=0.25$','$D_\Delta=0.5$','$D_\Delta=0.75$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, Watson_path, labels_system, figure_dir, ['/5_MeanDdelta_boxplot_Watson_fat_structures_Ddelta_SNR' num2str(SNR) '.pdf'], method_color, 'MeanDaniso')
% % 
% % 
%         %% Var[Ddelta] systems (Mean(D_delta^2) = 0.425)
%         system_name_list = ["bimodal_ddelta_0.6_0.7","bimodal_ddelta_0.5_0.77","bimodal_ddelta_0.4_0.83","bimodal_ddelta_0.3_0.87","bimodal_ddelta_0.2_0.9"];
%         labels_system = {'$D_\Delta=(0.6,0.7)$','$D_\Delta=(0.5,0.77)$','$D_\Delta=(0.4,0.83)$','$D_\Delta=(0.3,0.87)$','$D_\Delta=(0.2,0.9)$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, '', labels_system, figure_dir, ['/6_VarDdelta_boxplot_bimodal_ddelta_SNR' num2str(SNR) '.pdf'], method_color, 'VarDaniso')
%         
%         system_name_list = ["unimodal_variablevddeltan_0.01","unimodal_variablevddeltan_0.05","unimodal_variablevddeltan_0.1","unimodal_variablevddeltan_0.15"];
%         labels_system = {'V/m2=0.01','V/m2=0.05','V/m2=0.1','V/m2=0.15'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, '', labels_system, figure_dir, ['/6_VarDdelta_boxplot_unimodal_variablevddeltan_SNR' num2str(SNR) '.pdf'], method_color, 'VarDaniso')
% % 
% %        
% % 
% % %         %% Dispersive systems
%         system_name_list = ["Watson_fat_sticks_OP_0.01","Watson_fat_sticks_OP_0.25","Watson_fat_sticks_OP_0.5","Watson_fat_sticks_OP_0.75","single_fat_stick"];
%         labels_system = {'$\mathrm{OP}=0.01$','$\mathrm{OP}=0.25$','$\mathrm{OP}=0.5$','$\mathrm{OP}=0.75$','$\mathrm{OP}=1$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, Watson_path, labels_system, figure_dir, ['/7_Dispersion_boxplot_fat_SNR' num2str(SNR) '.pdf'], method_color)
%         
% 
%         system_name_list = ["Watson_fat_structures_OP_0.01","Watson_fat_structures_OP_0.25","Watson_fat_structures_OP_0.5","Watson_fat_structures_OP_0.75","unimodal_setvddeltan_0.6816"];
%         labels_system = {'$\mathrm{OP}=0.01$','$\mathrm{OP}=0.25$','$\mathrm{OP}=0.5$','$\mathrm{OP}=0.75$','$\mathrm{OP}=1$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, Watson_path, labels_system, figure_dir, ['/7_Dispersion_boxplot_structures_SNR' num2str(SNR) '.pdf'], method_color)
% % 
% 
% 
% 
% %         system_name_list = ["dispersive_real_stick_0.1","dispersive_real_stick_0.3","dispersive_real_stick_0.5","dispersive_real_stick_0.8","dispersive_real_stick_0.9"];
% %         labels_system = {'$\mathrm{OP}=0.67$','$\mathrm{OP}=0.8$','$\mathrm{OP}=0.9$','$\mathrm{OP}=0.96$','$\mathrm{OP}=0.98$'};
% %         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, Watson_path, labels_system, figure_dir, ['/7_Dispersion_boxplot_real_SNR' num2str(SNR) '.pdf'], method_color)
% %         
% %         % '0.1', OP=0.6736
% %         % '0.3', OP=0.7981
% %         % '0.5', OP=0.9018
% %         % '0.8', OP=0.9565
% %         % '0.9', OP=0.9837
% %         
% %         %% Crossing systems
%         system_name_list = ["two_real_sticks_crossing_15","two_real_sticks_crossing_30","two_real_sticks_crossing_45","two_real_sticks_crossing_70","two_real_sticks_crossing_90"];
%         labels_system = {'$\alpha = 15^\circ$','$\alpha = 30^\circ$','$\alpha = 45^\circ$','$\alpha = 70^\circ$','$\alpha = 90^\circ$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, Watson_path, labels_system, figure_dir, ['/8_Crossing_two-way_boxplot_real_SNR' num2str(SNR) '.pdf'], method_color)
% % 
%         system_name_list = ["three_real_sticks_crossing_15","three_real_sticks_crossing_30","three_real_sticks_crossing_45","three_real_sticks_crossing_70","three_real_sticks_crossing_90"];
%         labels_system = {'$\alpha = 15^\circ$','$\alpha = 30^\circ$','$\alpha = 45^\circ$','$\alpha = 70^\circ$','$\alpha = 90^\circ$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, Watson_path, labels_system, figure_dir, ['/8_Crossing_three-way_boxplot_real_SNR' num2str(SNR) '.pdf'], method_color)
% %         
        %% Partial volume [fibers,CSF] systems
        system_name_list = ["composite_real_stick_and_CSFf_0","composite_real_stick_and_CSFf_0.25","composite_real_stick_and_CSFf_0.5","composite_real_stick_and_CSFf_0.75","composite_real_stick_and_CSFf_1"];
        labels_system = {'$f_{\mathrm{iso}}=0$','$f_{\mathrm{iso}}=0.25$','$f_{\mathrm{iso}}=0.5$','$f_{\mathrm{iso}}=0.75$','$f_{\mathrm{iso}}=1$'};
        boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, Watson_path, labels_system, figure_dir, ['/9_Partial_volume_boxplot_real_stick_CSFf_SNR' num2str(SNR) '.pdf'], method_color)
% %         
% %         %% Partial volume [tumors,edema] systems
%         system_name_list = ["composite_edema_and_gliomaf_0","composite_edema_and_gliomaf_0.25","composite_edema_and_gliomaf_0.5","composite_edema_and_gliomaf_0.75","composite_edema_and_gliomaf_1"];
%         labels_system = {'$f_{\mathrm{g}}=0$','$f_{\mathrm{g}}=0.25$','$f_{\mathrm{g}}=0.5$','$f_{\mathrm{g}}=0.75$','$f_{\mathrm{g}}=1$'};
%         boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, Watson_path, labels_system, figure_dir, ['/10_Partial_volume_boxplot_edema_gliomaf_SNR' num2str(SNR) '.pdf'], method_color)
% %         
        total_duration = round(toc(timestart)/60,1)
    end
end