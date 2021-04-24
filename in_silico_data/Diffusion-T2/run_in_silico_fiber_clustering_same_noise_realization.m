% Run inversion pipeline

% Write using '...', not "..."

close all
clearvars

xps_file = fullfile(pwd, 'data_xps.mat');

run_spread_OP = 0;
run_crossing_angle = 0;
run_three_way_crossing = 0;
run_three_way_crossing_myelin = 1;

method = 'dtr2d';

%% Prepare simulations
if run_three_way_crossing
    structure_info.compartment_names = {'intra'; 'extra'};
    structure_info.method = method;
    structure_info.relative_std = [0.01 0.01];
    structure_info.fraction = [0.45, 0.55];
    structure_info.N = 50;
    
    inter_fiber_fraction = [1/3 1/3 1/3];
    
    case_name = 'intra_fiber_dist_three_way_crossing_more_iso';
    NBS_list = [100];
    do_add_iso = 1;
    
    config_n = 1;
    case_name = [case_name '_' num2str(config_n)]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADD HERE FOR OTHER NOISE REALIZATION %%%%%%%%%%%%%%%%%%
    
    if do_add_iso
        output_directory = fullfile(pwd, strcat(case_name, '_iso'));
    else
        output_directory = fullfile(pwd, strcat(case_name));
    end
    
    structure_info.output_dir = output_directory;
    
    if ~exist(output_directory, 'dir')
        msf_mkdir(output_directory);
    end
    
    % Each voxel contains a different numerical experiment
    SNR_list = 200*ones([1 8]); %[30 50 70 90 30 50 70 90];
%     SNR_list = 200*ones([1 8]);
    SNR_matrix = reshape(SNR_list, [2 2 2]);
    sz_SNR = size(SNR_matrix);
    
    theta_x_gt = pi/2;
    phi_x_gt = 0;
    
    theta_y_gt = pi/4; % Technically along y-z
    phi_y_gt = pi/2;
    
    theta_z_gt = 0;
    phi_z_gt = 0;
    
    struct_x = structure_info;
    struct_y = structure_info;
    struct_z = structure_info;
    struct_x.case_name = strcat('Config', num2str(config_n), '_fiber_x');
    struct_y.case_name = strcat('Config', num2str(config_n), '_fiber_yz');
    struct_z.case_name = strcat('Config', num2str(config_n), '_fiber_z');
    struct_x.mean_diso = [0.75, 1.1]*1e-9;
    struct_y.mean_diso = [0.75, 1]*1e-9;
    struct_z.mean_diso = [0.75, 0.9]*1e-9;
    struct_x.mean_ddelta = [0.95, 0.65];
    struct_y.mean_ddelta = [0.95, 0.7];
    struct_z.mean_ddelta = [0.95, 0.75];
    struct_x.mean_theta = theta_x_gt*[1 1];
    struct_y.mean_theta = theta_y_gt*[1 1];
    struct_z.mean_theta = theta_z_gt*[1 1];
    struct_x.mean_phi = phi_x_gt*[1 1];
    struct_y.mean_phi = phi_y_gt*[1 1];
    struct_z.mean_phi = phi_z_gt*[1 1];
    
    struct_x.mean_dpar = struct_x.mean_diso.*(1 + 2*struct_x.mean_ddelta);
    struct_y.mean_dpar = struct_y.mean_diso.*(1 + 2*struct_y.mean_ddelta);
    struct_z.mean_dpar = struct_z.mean_diso.*(1 + 2*struct_z.mean_ddelta);
    struct_x.mean_dperp = struct_x.mean_diso.*(1 - struct_x.mean_ddelta);
    struct_y.mean_dperp = struct_y.mean_diso.*(1 - struct_y.mean_ddelta);
    struct_z.mean_dperp = struct_z.mean_diso.*(1 - struct_z.mean_ddelta);
    
    if config_n == 1
        struct_x.mean_t = [75, 60]*1e-3;
        struct_y.mean_t = [100, 70]*1e-3;
        struct_z.mean_t = [120, 85]*1e-3;
        
    elseif config_n == 2
        struct_x.mean_t = [110, 85]*1e-3;
        struct_y.mean_t = [80, 65]*1e-3;
        struct_z.mean_t = [100, 70]*1e-3;
        
    elseif config_n == 3
        struct_x.mean_t = [95, 70]*1e-3;
        struct_y.mean_t = [120, 75]*1e-3;
        struct_z.mean_t = [70, 60]*1e-3;
        
%         struct_x.mean_t*[0.45 0.55]'*1e3
%         struct_y.mean_t*[0.45 0.55]'*1e3
%         struct_z.mean_t*[0.45 0.55]'*1e3
    end
    
    struct_x_out = create_heterogeneous_fiber(struct_x);
    struct_y_out = create_heterogeneous_fiber(struct_y);
    struct_z_out = create_heterogeneous_fiber(struct_z);
    
    dpar_x = struct_x_out.dpar;
    dpar_y = struct_y_out.dpar;
    dpar_z = struct_z_out.dpar;
    dperp_x = struct_x_out.dperp;
    dperp_y = struct_y_out.dperp;
    dperp_z = struct_z_out.dperp;
    theta_x = struct_x_out.theta;
    theta_y = struct_y_out.theta;
    theta_z = struct_z_out.theta;
    phi_x = struct_x_out.phi;
    phi_y = struct_y_out.phi;
    phi_z = struct_z_out.phi;
    r_x = struct_x_out.r;
    r_y = struct_y_out.r;
    r_z = struct_z_out.r;
    w_x = struct_x_out.w;
    w_y = struct_y_out.w;
    w_z = struct_z_out.w;
    
    diso_iso = 2*1e-9;
    t_iso = 500*1e-3;
    w_iso = 0.2;
    
    ground_truth = struct;
    ground_truth.dpar = [sum(dpar_x.*w_x) sum(dpar_y.*w_y) sum(dpar_z.*w_z)];
    ground_truth.dperp = [sum(dperp_x.*w_x) sum(dperp_y.*w_y) sum(dperp_z.*w_z)];
    ground_truth.theta = [theta_x_gt theta_y_gt theta_z_gt];
    ground_truth.phi = [phi_x_gt phi_y_gt phi_z_gt];
    ground_truth.r = [sum(r_x.*w_x) sum(r_y.*w_y) sum(r_z.*w_z)];
    ground_truth.w = inter_fiber_fraction;
    save(fullfile(output_directory, 'ground_truth_three_way_crossing.mat'), 'ground_truth');
    
    for n = 1:length(struct_x.mean_diso)
        ground_truth = struct;
        ground_truth.dpar = [struct_x.mean_dpar(n) struct_y.mean_dpar(n) struct_z.mean_dpar(n)];
        ground_truth.dperp = [struct_x.mean_dperp(n) struct_y.mean_dperp(n) struct_z.mean_dperp(n)];
        ground_truth.theta = [struct_x.mean_theta(n) struct_y.mean_theta(n) struct_z.mean_theta(n)];
        ground_truth.phi = [struct_x.mean_phi(n) struct_y.mean_phi(n) struct_z.mean_phi(n)];
        ground_truth.r = [1/struct_x.mean_t(n) 1/struct_y.mean_t(n) 1/struct_z.mean_t(n)];
        ground_truth.w = [struct_x.fraction(n) struct_y.fraction(n) struct_z.fraction(n)];
        save(fullfile(output_directory, strcat('ground_truth_three_way_crossing_', structure_info.compartment_names{n},'.mat')), 'ground_truth');
    end
    
    dpar = [dpar_x ; dpar_y ; dpar_z];
    dperp = [dperp_x ; dperp_y ; dperp_z];
    theta = [theta_x ; theta_y ; theta_z];
    phi = [phi_x ; phi_y ; phi_z];
    r = [r_x ; r_y ; r_z];
    w = inter_fiber_fraction.*[w_x ; w_y ; w_z];
    w = w/sum(w);
    
    load(xps_file, 'xps');
    fake_data = zeros([sz_SNR xps.n]);
    sz = size(fake_data);
    fake_mask = ones(sz(1:3));
    SNR_map = zeros(sz(1:3));

    for i_SNR_1 = 1:sz_SNR(1)
        for i_SNR_2 = 1:sz_SNR(2)
            for i_SNR_3 = 1:sz_SNR(3)
                
                vx = i_SNR_1;
                vy = i_SNR_2;
                vz = i_SNR_3;
                
                SNR = SNR_matrix(i_SNR_1, i_SNR_2, i_SNR_3);
                SNR_map(vx, vy, vz) = SNR;
                
                if do_add_iso
                    dpar_all = [dpar ; diso_iso];
                    dperp_all = [dperp ; diso_iso];
                    theta_all = [theta ; 0];
                    phi_all = [phi ; 0];
                    r_all = [r ; 1/t_iso];
                    w_all = [(1-w_iso)*w ; w_iso];
                    w_all = w_all/sum(w_all);
                else
                    dpar_all = dpar;
                    dperp_all = dperp;
                    theta_all = theta;
                    phi_all = phi;
                    r_all = r;
                    w_all = w;
                end
                
                if i_SNR_1 == 1 && i_SNR_2 == 1 && i_SNR_3 == 1
                    syst_struct.dpar = dpar_all;
                    syst_struct.dperp = dperp_all;
                    syst_struct.theta = theta_all;
                    syst_struct.phi = phi_all;
                    syst_struct.r = r_all;
                    syst_struct.w = w_all;
                    save(fullfile(output_directory, 'ground_truth_system.mat'), 'syst_struct');
                end
                
                % Simulate signal
                opt = mdm_opt();
                opt = dtr2d_opt(opt);
                opt.dtr2d.n_out = length(dpar_all); % Need to keep all the components !!!!!!!!!!!!!!!!!!!!!!!!!!!
                dtd = dtr2d_par2dist(dpar_all, dperp_all, theta_all, phi_all, r_all, w_all);
                m = dtr2d_dtr2d2m(dtd, opt);
                s_true = dtr2d_1d_fit2data(m, xps);
              
%                 b_theta = acos(xps.u(:,3)./sqrt(xps.u(:,1).^2 + xps.u(:,2).^2 + xps.u(:,3).^2));
%                 b_phi = atan2(xps.u(:,2),xps.u(:,1));
%                 if strcmp(method, 'dtd')
%                     acq_mtrx = [round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
%                 elseif strcmp(method, 'dtr1d')
%                     acq_mtrx = [xps.tr round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
%                 elseif strcmp(method, 'dtr2d')
%                     acq_mtrx = [xps.te round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
%                 elseif strcmp(method, 'dtr2r1d')
%                     acq_mtrx = [xps.tr xps.te round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
%                 end
%                 [acq_mtrx_sort,indx_sort] = sortrows(acq_mtrx);
                
%                 figure()
%                 plot(1:xps.n, s_true(indx_sort))
                
%                 list_TE = unique(xps.te);
%                 N_TE = length(list_TE);
%                 s = zeros([xps.n 1]);
%                 for i = 1:N_TE
%                    ind = xps.te == list_TE(i);
%                    i_b0 = find(xps.te == list_TE(i) & xps.b == 0, 1, 'first');
%                    max_s = s_true(i_b0);
%                    
%                    if isfinite(SNR)
%                        s = sqrt((s_true + 1/SNR*randn([xps.n 1])).^2 + (1/SNR*randn([xps.n 1])).^2); % Adding Rician noise
%                    else
%                        s = s_true;
%                    end
%                 end
                
                if isfinite(SNR)
                    s = sqrt((s_true + 1/SNR*randn([xps.n 1])).^2 + (1/SNR*randn([xps.n 1])).^2); % Adding Rician noise
                else
                    s = s_true;
                end
                
                fake_data(vx, vy, vz, :) = s;
            end
        end
    end
    
    mdm_nii_write(SNR_map, fullfile(output_directory, 'SNR_map.nii.gz'));
    
elseif run_three_way_crossing_myelin
    f_myelin = 0.1111;
    
    structure_info.compartment_names = {'intra'; 'extra'; 'myelin'};
    structure_info.method = method;
    structure_info.relative_std = [0.01 0.01 0.01];
    structure_info.fraction = [0.45*(1-f_myelin), 0.55*(1-f_myelin), f_myelin];
    structure_info.N = 50;
    
    inter_fiber_fraction = [1/3 1/3 1/3];
    
    case_name = 'intra_fiber_dist_three_way_crossing_myelin_GT';
    NBS_list = [100];
    do_add_iso = 1;
    
    config_n = 1;
    case_name = [case_name '_' num2str(config_n)]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADD HERE FOR OTHER NOISE REALIZATION %%%%%%%%%%%%%%%%%%
    
    if do_add_iso
        output_directory = fullfile(pwd, strcat(case_name, '_iso'));
    else
        output_directory = fullfile(pwd, strcat(case_name));
    end
    
    structure_info.output_dir = output_directory;
    
    if ~exist(output_directory, 'dir')
        msf_mkdir(output_directory);
    end
    
    % Each voxel contains a different numerical experiment
    SNR_list = 200*ones([1 8]); %[30 50 70 90 30 50 70 90]; 200*ones([1 8]);
%     SNR_list = 200*ones([1 8]);
    SNR_matrix = reshape(SNR_list, [2 2 2]);
    sz_SNR = size(SNR_matrix);
    
    theta_x_gt = pi/2;
    phi_x_gt = 0;
    
    theta_y_gt = pi/4; % Technically along y-z
    phi_y_gt = pi/2;
    
    theta_z_gt = 0;
    phi_z_gt = 0;
    
    struct_x = structure_info;
    struct_y = structure_info;
    struct_z = structure_info;
    struct_x.case_name = strcat('Config', num2str(config_n), '_fiber_x');
    struct_y.case_name = strcat('Config', num2str(config_n), '_fiber_yz');
    struct_z.case_name = strcat('Config', num2str(config_n), '_fiber_z');
    struct_x.mean_diso = [0.75, 1.1, 0.4]*1e-9;
    struct_y.mean_diso = [0.75, 1, 0.4]*1e-9;
    struct_z.mean_diso = [0.75, 0.9, 0.4]*1e-9;
    struct_x.mean_ddelta = [0.95, 0.65, 0.95];
    struct_y.mean_ddelta = [0.95, 0.7, 0.95];
    struct_z.mean_ddelta = [0.95, 0.75, 0.95];
    struct_x.mean_theta = theta_x_gt*[1 1 1];
    struct_y.mean_theta = theta_y_gt*[1 1 1];
    struct_z.mean_theta = theta_z_gt*[1 1 1];
    struct_x.mean_phi = phi_x_gt*[1 1 1];
    struct_y.mean_phi = phi_y_gt*[1 1 1];
    struct_z.mean_phi = phi_z_gt*[1 1 1];
    
    struct_x.mean_dpar = struct_x.mean_diso.*(1 + 2*struct_x.mean_ddelta);
    struct_y.mean_dpar = struct_y.mean_diso.*(1 + 2*struct_y.mean_ddelta);
    struct_z.mean_dpar = struct_z.mean_diso.*(1 + 2*struct_z.mean_ddelta);
    struct_x.mean_dperp = struct_x.mean_diso.*(1 - struct_x.mean_ddelta);
    struct_y.mean_dperp = struct_y.mean_diso.*(1 - struct_y.mean_ddelta);
    struct_z.mean_dperp = struct_z.mean_diso.*(1 - struct_z.mean_ddelta);
    
    if config_n == 1
        struct_x.mean_t = [75, 60, 10]*1e-3;
        struct_y.mean_t = [100, 70, 10]*1e-3;
        struct_z.mean_t = [120, 85, 10]*1e-3;
        
    elseif config_n == 2
        struct_x.mean_t = [110, 85, 10]*1e-3;
        struct_y.mean_t = [80, 65, 10]*1e-3;
        struct_z.mean_t = [100, 70, 10]*1e-3;
        
    elseif config_n == 3
        struct_x.mean_t = [95, 70, 10]*1e-3;
        struct_y.mean_t = [120, 75, 10]*1e-3;
        struct_z.mean_t = [70, 60, 10]*1e-3;
        
%         struct_x.mean_t*[0.45 0.55]'*1e3
%         struct_y.mean_t*[0.45 0.55]'*1e3
%         struct_z.mean_t*[0.45 0.55]'*1e3
    end
    
    struct_x_out = create_heterogeneous_fiber(struct_x);
    struct_y_out = create_heterogeneous_fiber(struct_y);
    struct_z_out = create_heterogeneous_fiber(struct_z);
    
    dpar_x = struct_x_out.dpar;
    dpar_y = struct_y_out.dpar;
    dpar_z = struct_z_out.dpar;
    dperp_x = struct_x_out.dperp;
    dperp_y = struct_y_out.dperp;
    dperp_z = struct_z_out.dperp;
    theta_x = struct_x_out.theta;
    theta_y = struct_y_out.theta;
    theta_z = struct_z_out.theta;
    phi_x = struct_x_out.phi;
    phi_y = struct_y_out.phi;
    phi_z = struct_z_out.phi;
    r_x = struct_x_out.r;
    r_y = struct_y_out.r;
    r_z = struct_z_out.r;
    w_x = struct_x_out.w;
    w_y = struct_y_out.w;
    w_z = struct_z_out.w;
    
    diso_iso = 2*1e-9;
    t_iso = 500*1e-3;
    w_iso = 0.1;
    
    ground_truth = struct;
    ground_truth.dpar = [sum(dpar_x.*w_x) sum(dpar_y.*w_y) sum(dpar_z.*w_z)];
    ground_truth.dperp = [sum(dperp_x.*w_x) sum(dperp_y.*w_y) sum(dperp_z.*w_z)];
    ground_truth.theta = [theta_x_gt theta_y_gt theta_z_gt];
    ground_truth.phi = [phi_x_gt phi_y_gt phi_z_gt];
    ground_truth.r = [sum(r_x.*w_x) sum(r_y.*w_y) sum(r_z.*w_z)];
    ground_truth.w = inter_fiber_fraction;
    save(fullfile(output_directory, 'ground_truth_three_way_crossing.mat'), 'ground_truth');
    
    for n = 1:length(struct_x.mean_diso)
        ground_truth = struct;
        ground_truth.dpar = [struct_x.mean_dpar(n) struct_y.mean_dpar(n) struct_z.mean_dpar(n)];
        ground_truth.dperp = [struct_x.mean_dperp(n) struct_y.mean_dperp(n) struct_z.mean_dperp(n)];
        ground_truth.theta = [struct_x.mean_theta(n) struct_y.mean_theta(n) struct_z.mean_theta(n)];
        ground_truth.phi = [struct_x.mean_phi(n) struct_y.mean_phi(n) struct_z.mean_phi(n)];
        ground_truth.r = [1/struct_x.mean_t(n) 1/struct_y.mean_t(n) 1/struct_z.mean_t(n)];
        ground_truth.w = [struct_x.fraction(n) struct_y.fraction(n) struct_z.fraction(n)];
        save(fullfile(output_directory, strcat('ground_truth_three_way_crossing_', structure_info.compartment_names{n},'.mat')), 'ground_truth');
    end
    
    dpar = [dpar_x ; dpar_y ; dpar_z];
    dperp = [dperp_x ; dperp_y ; dperp_z];
    theta = [theta_x ; theta_y ; theta_z];
    phi = [phi_x ; phi_y ; phi_z];
    r = [r_x ; r_y ; r_z];
    w = inter_fiber_fraction.*[w_x ; w_y ; w_z];
    w = w/sum(w);
    
    load(xps_file, 'xps');
    fake_data = zeros([sz_SNR xps.n]);
    sz = size(fake_data);
    fake_mask = ones(sz(1:3));
    SNR_map = zeros(sz(1:3));

    for i_SNR_1 = 1:sz_SNR(1)
        for i_SNR_2 = 1:sz_SNR(2)
            for i_SNR_3 = 1:sz_SNR(3)
                
                vx = i_SNR_1;
                vy = i_SNR_2;
                vz = i_SNR_3;
                
                SNR = SNR_matrix(i_SNR_1, i_SNR_2, i_SNR_3);
                SNR_map(vx, vy, vz) = SNR;
                
                if do_add_iso
                    dpar_all = [dpar ; diso_iso];
                    dperp_all = [dperp ; diso_iso];
                    theta_all = [theta ; 0];
                    phi_all = [phi ; 0];
                    r_all = [r ; 1/t_iso];
                    w_all = [(1-w_iso)*w ; w_iso];
                    w_all = w_all/sum(w_all);
                else
                    dpar_all = dpar;
                    dperp_all = dperp;
                    theta_all = theta;
                    phi_all = phi;
                    r_all = r;
                    w_all = w;
                end
                
                if i_SNR_1 == 1 && i_SNR_2 == 1 && i_SNR_3 == 1
                    syst_struct.dpar = dpar_all;
                    syst_struct.dperp = dperp_all;
                    syst_struct.theta = theta_all;
                    syst_struct.phi = phi_all;
                    syst_struct.r = r_all;
                    syst_struct.w = w_all;
                    save(fullfile(output_directory, 'ground_truth_system.mat'), 'syst_struct');
                end
                
                % Simulate signal
                opt = mdm_opt();
                opt = dtr2d_opt(opt);
                opt.dtr2d.n_out = length(dpar_all); % Need to keep all the components !!!!!!!!!!!!!!!!!!!!!!!!!!!
                dtd = dtr2d_par2dist(dpar_all, dperp_all, theta_all, phi_all, r_all, w_all);
                m = dtr2d_dtr2d2m(dtd, opt);
                s_true = dtr2d_1d_fit2data(m, xps);
              
%                 b_theta = acos(xps.u(:,3)./sqrt(xps.u(:,1).^2 + xps.u(:,2).^2 + xps.u(:,3).^2));
%                 b_phi = atan2(xps.u(:,2),xps.u(:,1));
%                 if strcmp(method, 'dtd')
%                     acq_mtrx = [round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
%                 elseif strcmp(method, 'dtr1d')
%                     acq_mtrx = [xps.tr round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
%                 elseif strcmp(method, 'dtr2d')
%                     acq_mtrx = [xps.te round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
%                 elseif strcmp(method, 'dtr2r1d')
%                     acq_mtrx = [xps.tr xps.te round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
%                 end
%                 [acq_mtrx_sort,indx_sort] = sortrows(acq_mtrx);
                
%                 figure()
%                 plot(1:xps.n, s_true(indx_sort))
                
%                 list_TE = unique(xps.te);
%                 N_TE = length(list_TE);
%                 s = zeros([xps.n 1]);
%                 for i = 1:N_TE
%                    ind = xps.te == list_TE(i);
%                    i_b0 = find(xps.te == list_TE(i) & xps.b == 0, 1, 'first');
%                    max_s = s_true(i_b0);
%                    
%                    if isfinite(SNR)
%                        s = sqrt((s_true + 1/SNR*randn([xps.n 1])).^2 + (1/SNR*randn([xps.n 1])).^2); % Adding Rician noise
%                    else
%                        s = s_true;
%                    end
%                 end
                
                if isfinite(SNR)
                    s = sqrt((s_true + 1/SNR*randn([xps.n 1])).^2 + (1/SNR*randn([xps.n 1])).^2); % Adding Rician noise
                else
                    s = s_true;
                end
                
                fake_data(vx, vy, vz, :) = s;
            end
        end
    end
    
    mdm_nii_write(SNR_map, fullfile(output_directory, 'SNR_map.nii.gz'));
end

for NBS = NBS_list
    fprintf(strcat('####################################################################################### Nbs=', num2str(NBS), '\n'))
    
    output_directory_Nbs = fullfile(output_directory, strcat('Nbs', num2str(NBS)));
    
    if ~exist(output_directory_Nbs, 'dir')
        msf_mkdir(output_directory_Nbs);
    end

    mdm_nii_write(fake_data, fullfile(output_directory_Nbs, 'data.nii.gz'));
    [~, h] = mdm_nii_read(fullfile(output_directory_Nbs, 'data.nii.gz'));
    mdm_nii_write(fake_mask, fullfile(output_directory_Nbs, 'data_mask.nii.gz'), h);
    copyfile(xps_file, output_directory_Nbs);
    
    input_parameters = struct;
    
    %% Input files
    input_parameters.initial_directory = output_directory_Nbs; % Where the data can be found
    input_parameters.data_file = 'data.nii.gz'; % Diffusion/relaxation data
    input_parameters.mask_file = 'data_mask.nii.gz'; % Data mask - If no mask, simply write '', one will be generated from mdm_s_mask
    input_parameters.xps_file = 'data_xps.mat'; % xps matfile
    
    % input_parameters.initial_directory = fullfile(pwd,'Data_test'); % Where the data can be found
    % input_parameters.mask_file = 'b0_bet_mask_data.nii.gz'; % Data mask - If no mask, simply write '', one will be generated from mdm_s_mask
    
    %% Inversion method
    input_parameters.inversion_method = 'dtr2d';
    input_parameters.do_covariance = 0;

    %% Data correction (extrapolated references + Elastix)
    input_parameters.do_data_correction = 0;
    input_parameters.do_rotate_bvec = 0;
    
    %% Global parameters
    input_parameters.nb_MC_inversions = NBS; % Number of Monte Carlo realizations
    input_parameters.nb_parallel_workers = Inf; % Number of parallel workers, put Inf to have the largest number of workers on your machine
    input_parameters.threshold_d_iso_big = 2.5; % Threshold on d_iso to define the big bin (in {\mu}m^2/ms)
    input_parameters.threshold_d_delta_thin = 0.5; % Threshold on anisotropy to define the thin bin
    input_parameters.threshold_weight_thin = 0.05; % Threshold on the thin-bin total weight to create a thin-bin data mask (useless for now)
    input_parameters.dir_flips = [0 0 0];
    
    %% Colormap bounds
    input_parameters.clim.diso_clim = 3.5e-9*[0 1];
    input_parameters.clim.msqddelta_clim = 1*[0 1];
    input_parameters.clim.vdiso_clim = .3*3e-9^2*[0 1];
    input_parameters.clim.vsqddelta_clim = .1*[0 1];
    input_parameters.clim.cvdisosqddelta_clim = .1*3e-9*1*[-1 1];
    
    if strcmp(input_parameters.inversion_method, 'dtr2d')
        mr_clim = 30*[0 1];
        input_parameters.clim.mr_clim = mr_clim;
        input_parameters.clim.vr_clim = .2*max(mr_clim)^2*[0 1];
        input_parameters.clim.cvdisor_clim = 0.1*max(mr_clim)*3e-9*[-1 1];
        input_parameters.clim.cvsqddeltar_clim = 0.1*max(mr_clim)*[-1 1];
        
        input_parameters.clim.mt_clim = [1/max(mr_clim) 0.11];
        input_parameters.clim.mt_clim_blackandwhite = [1/max(mr_clim) 0.15];
        input_parameters.clim.vt_clim = 0.5*[0 1];
        input_parameters.clim.cvdisot_clim = 1e-9*[-1 1];
        input_parameters.clim.cvsqddeltat_clim = 0.05*[-1 1];
        
    elseif strcmp(input_parameters.inversion_method, 'dtr1d')
        mr_clim = 0.8*[0 1];
        input_parameters.clim.mr_clim = mr_clim;
        input_parameters.clim.vr_clim = .2*max(mr_clim)^2*[0 1];
        input_parameters.clim.cvdisor_clim = 0.1*max(mr_clim)*3e-9*[-1 1];
        input_parameters.clim.cvsqddeltar_clim = 0.1*max(mr_clim)*1*[-1 1];
        
        mt_clim = [0.5 5];
        input_parameters.clim.mt_clim = mt_clim;
        input_parameters.clim.mt_clim_blackandwhite = [1 10];
        input_parameters.clim.vt_clim = 0.25*[0 max(mt_clim)^2];
        input_parameters.clim.cvdisot_clim = 0.1*max(mt_clim)*3e-9*[-1 1];
        input_parameters.clim.cvsqddeltat_clim = 0.1*max(mt_clim)*1*[-1 1];
        
    elseif strcmp(input_parameters.inversion_method, 'dtr2r1d')
        mr2_clim = 30*[0 1];
        input_parameters.clim.mr2_clim = mr2_clim;
        input_parameters.clim.vr2_clim = .2*max(mr2_clim)^2*[0 1];
        input_parameters.clim.cvdisor2_clim = 0.1*max(mr2_clim)*3e-9*[-1 1];
        input_parameters.clim.cvsqddeltar2_clim = 0.1*max(mr2_clim)*[-1 1];
        
        mt2_clim = [1/max(mr2_clim) 0.11];
        input_parameters.clim.mt2_clim = mt2_clim;
        input_parameters.clim.mt2_clim_blackandwhite = [1/max(mr2_clim) 0.15];
        input_parameters.clim.vt2_clim = 0.5*[0 1];
        input_parameters.clim.cvdisot2_clim = 1e-9*[-1 1];
        input_parameters.clim.cvsqddeltat2_clim = 0.05*[-1 1];
        
        mr1_clim = 0.8*[0 1];
        input_parameters.clim.mr1_clim = mr1_clim;
        input_parameters.clim.vr1_clim = .2*max(mr1_clim)^2*[0 1];
        input_parameters.clim.cvdisor1_clim = 0.1*max(mr1_clim)*3e-9*[-1 1];
        input_parameters.clim.cvsqddeltar1_clim = 0.1*max(mr1_clim)*1*[-1 1];
        
        mt1_clim = [0.5 5];
        input_parameters.clim.mt1_clim = mt1_clim;
        input_parameters.clim.mt1_clim_blackandwhite = [1 10];
        input_parameters.clim.vt1_clim = 0.25*[0 max(mt1_clim)^2];
        input_parameters.clim.cvdisot1_clim = 0.1*max(mt1_clim)*3e-9*[-1 1];
        input_parameters.clim.cvsqddeltat1_clim = 0.1*max(mt1_clim)*1*[-1 1];
        
        input_parameters.clim.cvr1r2_clim = 0.1*max(mr1_clim)*max(mr2_clim)*[-1 1];
        input_parameters.clim.cvt1t2_clim = 0.1*max(mt1_clim)*max(mt2_clim)*[-1 1];
    end
    
    % If one wants to optimize the colormap bounds without deleting input_parameters.parameter_maps_directory everytime
    input_parameters.repeat_colormaps = 0;
    
    %% Do ODFs (and clustering)?
    input_parameters.do_ODFs = 1; % Do you want ODFs to be produced?
    input_parameters.framework_directory = '/Users/alexis_reymbaut/Dropbox/Research_Lund/Pipeline_general/Pipeline/md-dmri'; % Directory wherein 'methods' and 'tools' are located
    input_parameters.mrtrix_directory = ''; % Directory wherein the binaries of mrtrix are located
    input_parameters.max_nb_odf_peaks = 4; % Maximal number of ODF peaks per voxel
    input_parameters.threshold_ODF_w = 0.1; % 0.05, 0.1, 0.15, take out the core of the ODF for peak calculation
    input_parameters.nb_mesh_nodes = 1000; % 250, 350, 500, 1000, 3994, or 15970, to serve as projecting mesh for the discrete ODFs
    
    input_parameters.do_clustering = 1; % Do you want clustering to be performed? (Requires input_parameters.do_ODFs = 1 too.)
    struct_clustering.max_nb_clusters = input_parameters.max_nb_odf_peaks; % Maximal number of clusters
    struct_clustering.do_plot = 0; % Print clustering plots, do not use for big volumes
    struct_clustering.min_weight_per_cluster = 0.1; % Repeat clustering at the bootstrap level if one cluster does not weigh at least this
    input_parameters.structure_clustering = struct_clustering;
    
    %% Output directories
    input_parameters.data_directory = fullfile(input_parameters.initial_directory, '0_data');
    input_parameters.bootstrap_directory = fullfile(input_parameters.initial_directory, '1_bootstrap');
    input_parameters.parameter_maps_directory = fullfile(input_parameters.initial_directory, '2_nii_pdf_maps');
    input_parameters.odf_directory = fullfile(input_parameters.initial_directory, '3_odfs');
    input_parameters.clustering_directory = fullfile(input_parameters.initial_directory, '4_clustering');
    
    %% Run the inversion pipeline
    global_inversion_pipeline(input_parameters)
end