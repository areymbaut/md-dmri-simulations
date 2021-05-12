function [dpar, dperp, theta, phi, w] = choose_simulated_system(system_name, path, is_plot)

if nargin < 3
    is_plot=1;
end

if nargin < 2
    path='';
end

% "Fat sticks" are prolate tensors characterized by Diso = 0.8 and FA = 0.85
if contains(system_name,'fat')
    [dpar_fat, dperp_fat] = Dpara_Dperp_from_Diso_FA(0.8,0.85,'prolate');
end

% "Real sticks" follow Gaussian distributions from the get_real_stick function
    function [dpar, dperp, theta, phi, w] = set_of_real_sticks(mean_theta, mean_phi, path)
        N = length(mean_theta);
        dpar = [];
        dperp = [];
        theta = [];
        phi = [];
        w = [];
        for it = 1:N
            [dpar_list, dperp_list, theta_list, phi_list, w_list] = get_real_stick(mean_theta(it), mean_phi(it), path);
            dpar = [dpar ; dpar_list];
            dperp = [dperp ; dperp_list];
            theta = [theta ; theta_list];
            phi = [phi ; phi_list];
            w = [w ; w_list];
        end
        w = w/sum(w);
    end

% Microscopically isotropic systems

if startsWith(system_name,'unimodal_setvddeltan')
    vddeltan_value = 0.01;
    Diso = 0.8*1e-9;
    mean_value = str2double(extractAfter(system_name,'unimodal_setvddeltan_'));
    sigma_value = sqrt(vddeltan_value)*mean_value;
    ddelta = linspace(-0.5,1,500)';
    dpar = Diso*(1+2*ddelta);
    dperp = Diso*(1-ddelta);
    theta = zeros(size(dpar));
    phi = theta;
    w = exp(-1/2*(ddelta-mean_value).^2/sigma_value^2); w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(ddelta, w);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '.pdf'))
        clf(f)
    end
    
elseif startsWith(system_name,'unimodal_variablevddeltan')
    vddeltan_value = str2double(extractAfter(system_name,'unimodal_variablevddeltan_'));
    Diso = 0.8*1e-9;
    mean_value = 0.6519;
    sigma_value = sqrt(vddeltan_value)*mean_value;
    ddelta = linspace(-0.5,1,500)';
    dpar = Diso*(1+2*ddelta);
    dperp = Diso*(1-ddelta);
    theta = zeros(size(dpar));
    phi = theta;
    w = exp(-1/2*(ddelta-mean_value).^2/sigma_value^2); w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(ddelta, w);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '.pdf'))
        clf(f)
    end
    
    
elseif startsWith(system_name,'unimodal_setvdison')
    vdison_value = 0.01;
    D1 = str2double(extractAfter(system_name,'unimodal_setvdison_'));
    mean_value = D1*1e-9;
    sigma_value = sqrt(vdison_value)*mean_value;
    width_value = 3*sigma_value;
    dpar = linspace(max([0.01*1e-9,mean_value-width_value]),min([3*1e-9,mean_value+width_value]),200)';
    dperp = dpar;
    theta = zeros(size(dpar));
    phi = theta;
    w = exp(-1/2*(dpar-mean_value).^2/sigma_value^2); w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(dpar, w);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '.pdf'))
        clf(f)
    end
    
elseif startsWith(system_name,'unimodal_variablevdison')
    vdison_value = str2double(extractAfter(system_name,'unimodal_variablevdison_'));
    mean_value = 0.8*1e-9;
    sigma_value = sqrt(vdison_value)*mean_value;
    width_value = 3*sigma_value;
    dpar = linspace(max([0.01*1e-9,mean_value-width_value]),min([3*1e-9,mean_value+width_value]),200)';
    dperp = dpar;
    theta = zeros(size(dpar));
    phi = theta;
    w = exp(-1/2*(dpar-mean_value).^2/sigma_value^2); w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(dpar, w);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '.pdf'))
        clf(f)
    end
 
elseif startsWith(system_name,'bimodal_ddelta') 
    D1_char = extractBetween(system_name,'bimodal_ddelta_','_');
    D1_char = D1_char{1};
    D1 = str2double(D1_char);
    D2 = str2double(extractAfter(system_name,['bimodal_ddelta_' D1_char '_']));
    mean_values = [D1; D2];
    mean_value = sort(mean_values);
    sigma_value = 0.01;
    ddelta = linspace(-0.5,1,500)';
    diso = 0.8*1e-9;
    dpar = diso*(1+2*ddelta);
    dperp = diso*(1-ddelta);
    theta = zeros(size(dpar));
    phi = theta;
    w = exp(-1/2*(ddelta-mean_value(1)).^2/sigma_value^2)+exp(-1/2*(ddelta-mean_value(2)).^2/sigma_value^2); w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(ddelta, w);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '.pdf'))
        clf(f)
    end
    
elseif startsWith(system_name,'bimodal')
    D1_char = extractBetween(system_name,'bimodal_','_');
    D1_char = D1_char{1};
    D1 = str2double(D1_char);
    D2 = str2double(extractAfter(system_name,['bimodal_' D1_char '_']));
    mean_values = [D1; D2]*1e-9;
    mean_value = sort(mean_values);
    sigma_value = 0.05*1e-9;
    width_value = 0.4*1e-9;
    dpar = linspace(max([0.01*1e-9,mean_value(1)-width_value]),min([3*1e-9,mean_value(2)+width_value]),500)';
    dperp = dpar;
    theta = zeros(size(dpar));
    phi = theta;
    w = exp(-1/2*(dpar-mean_value(1)).^2/sigma_value^2)+exp(-1/2*(dpar-mean_value(2)).^2/sigma_value^2); w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(dpar, w);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '.pdf'))
        clf(f)
    end
    

    
    
% Microscopically anisotropic systems

elseif strcmp(system_name,'single_fat_stick')
    dpar = dpar_fat*1e-9;
    dperp = dperp_fat*1e-9;
    theta = zeros(size(dpar));
    phi = theta;
    w = 1; w = w/sum(w);
    
elseif strcmp(system_name,'single_real_stick')
    [dpar, dperp, theta, phi, w] = get_real_stick(0, 0, path);
    
elseif startsWith(system_name,'dispersive_real_stick')
    OP_char = extractAfter(system_name,'dispersive_real_stick_');
    [dpar, dperp, theta, phi, w] = get_real_stick_dispersive(0, 0, OP_char, path);
    %
    if is_plot
        [dpar_sorted, ind_dpar_sorted] = sort(dpar);
        w_dpar_sorted = w(ind_dpar_sorted);
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(dpar_sorted, w_dpar_sorted);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '_dpar.pdf'))
        clf(f)
        %
        [dperp_sorted, ind_dperp_sorted] = sort(dperp);
        w_dperp_sorted = w(ind_dperp_sorted);
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(dperp_sorted, w_dperp_sorted);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '_dperp.pdf'))
        clf(f)
        %
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        cla
        ax = polaraxes;
        for i = 1:length(phi)
            polarscatter(ax, phi(i), sin(theta(i)), 'b', 'filled', 'MarkerFaceAlpha', w(i)/max(w))
            hold on
            polarscatter(ax, phi(i), sin(theta(i)), 'r', 'filled', 'MarkerFaceAlpha', (1-w(i)/max(w))/10)
            hold on
        end
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'top';
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '_angular.pdf'))
        cla
        clf(f)
    end
        
elseif strcmp(system_name,'fat_stick_powder')
    Ncomp = 100; 
    dpar = dpar_fat*1e-9*ones([Ncomp 1]);
    dperp = dperp_fat*1e-9*ones([Ncomp 1]);
    load(fullfile(path,num2str(Ncomp)),'theta','phi');
    w = ones([Ncomp 1]); w = w/sum(w);
    
elseif strcmp(system_name,'real_stick_powder')
    Ncomp = 100; 
    load(fullfile(path,num2str(Ncomp)),'theta','phi');
    mean_theta = theta;
    mean_phi = phi;
    [dpar, dperp, theta, phi, w] = set_of_real_sticks(mean_theta, mean_phi, path);
    
% Crossing systems

elseif startsWith(system_name,'two_fat_sticks_crossing')
    crossing_angle_char = extractAfter(system_name,'two_fat_sticks_crossing_');
    crossing_angle = deg2rad(str2double(crossing_angle_char));
    dpar = dpar_fat*1e-9*ones([2 1]);
    dperp = dperp_fat*1e-9*ones([2 1]);
    theta = [pi/2 pi/2]';
    phi = [0 crossing_angle]';
    w = [1 1]'; w = w/sum(w);
    
elseif startsWith(system_name,'two_real_sticks_crossing')
    crossing_angle_char = extractAfter(system_name,'two_real_sticks_crossing_');
    crossing_angle = deg2rad(str2double(crossing_angle_char));
    mean_theta = [pi/2 pi/2]';
    mean_phi = [0 crossing_angle]';
    [dpar, dperp, theta, phi, w] = set_of_real_sticks(mean_theta, mean_phi, path);
   
   
elseif startsWith(system_name,'three_fat_sticks_crossing')
    crossing_angle_char = extractAfter(system_name,'three_fat_sticks_crossing_');
    crossing_angle = deg2rad(str2double(crossing_angle_char));
    dpar = dpar_fat*1e-9*ones([3 1]);
    dperp = dperp_fat*1e-9*ones([3 1]);
    theta = [pi/2 pi/2 asin(cos(crossing_angle)/cos(crossing_angle/2))]';
    phi = [0 crossing_angle crossing_angle/2]';
    w = [1 1 1]'; w = w/sum(w);

elseif startsWith(system_name,'three_real_sticks_crossing')
    crossing_angle_char = extractAfter(system_name,'three_real_sticks_crossing_');
    crossing_angle = deg2rad(str2double(crossing_angle_char));
    mean_theta = [pi/2 pi/2 asin(cos(crossing_angle)/cos(crossing_angle/2))]';
    mean_phi = [0 crossing_angle crossing_angle/2]';
    [dpar, dperp, theta, phi, w] = set_of_real_sticks(mean_theta, mean_phi, path);
    

    
    
% Dispersive systems
    
elseif startsWith(system_name,'Watson_fat_structures_Ddelta')
    vddeltan_value = 0.01;
    Diso = 0.8*1e-9;
    mean_value = str2double(extractAfter(system_name,'Watson_fat_structures_Ddelta_'));
    sigma_value = sqrt(vddeltan_value)*mean_value;
    ddelta = linspace(-0.5,1,500)';
    dpar = Diso*(1+2*ddelta);
    dperp = Diso*(1-ddelta);
    w = exp(-1/2*(ddelta-mean_value).^2/sigma_value^2);
    load(fullfile(path,'theta_phi_Watson_OP_0.01.mat'),'theta','phi');
    phi = phi'; % TO CORRECT
    
    Ncomp = length(dpar);
    dpar_all = [];
    dperp_all = [];
    theta_all = [];
    phi_all = [];
    w_all = [];
    for i = 1:length(theta) 
        dpar_all = [dpar_all ; dpar];
        dperp_all = [dperp_all ; dperp];
        theta_all = [theta_all ; theta(i)*ones([Ncomp 1])];
        phi_all = [phi_all ; phi(i)*ones([Ncomp 1])];
        w_all = [w_all ; w];
    end
    
    dpar = dpar_all;
    dperp = dperp_all;
    theta = theta_all;
    phi = phi_all;
    w = w_all/sum(w_all);
    
elseif startsWith(system_name,'Watson_fat_structures_OP')
    vddeltan_value = 0.01;
    Diso = 0.8*1e-9;
    mean_value = 0.6816; % FA = 0.85
    sigma_value = sqrt(vddeltan_value)*mean_value;
    ddelta = linspace(-0.5,1,500)';
    dpar = Diso*(1+2*ddelta);
    dperp = Diso*(1-ddelta);
    w = exp(-1/2*(ddelta-mean_value).^2/sigma_value^2);
    OP = extractAfter(system_name,'Watson_fat_structures_OP_');
    load(fullfile(path,strcat('theta_phi_Watson_OP_', OP, '.mat')),'theta','phi');
    phi = phi'; % TO CORRECT
    
    Ncomp = length(dpar);
    dpar_all = [];
    dperp_all = [];
    theta_all = [];
    phi_all = [];
    w_all = [];
    for i = 1:length(theta) 
        dpar_all = [dpar_all ; dpar];
        dperp_all = [dperp_all ; dperp];
        theta_all = [theta_all ; theta(i)*ones([Ncomp 1])];
        phi_all = [phi_all ; phi(i)*ones([Ncomp 1])];
        w_all = [w_all ; w];
    end
    
    dpar = dpar_all;
    dperp = dperp_all;
    theta = theta_all;
    phi = phi_all;
    w = w_all/sum(w_all);
    
% elseif startsWith(system_name,'Watson_fat_structures_Ddelta')
%     D_Delta = str2double(extractAfter(system_name,'Watson_fat_structures_Ddelta_'));
%     load(fullfile(path,'theta_phi_Watson_OP_0.01.mat'),'theta','phi');
%     phi = phi'; % TO CORRECT
%     Ncomp = length(theta);
%     Diso = 0.8*1e-9;
%     dpar = Diso*(1+2*D_Delta)*ones([Ncomp 1]);
%     dperp = Diso*(1-D_Delta)*ones([Ncomp 1]);
%     w = ones([Ncomp 1]); w = w/sum(w);    
    
elseif startsWith(system_name,'Watson_fat_sticks_OP')
    OP = extractAfter(system_name,'Watson_fat_sticks_OP_');
    load(fullfile(path,strcat('theta_phi_Watson_OP_', OP, '.mat')),'theta','phi');
    phi = phi'; % TO CORRECT
    Ncomp = length(theta);
    dpar = dpar_fat*1e-9*ones([Ncomp 1]);
    dperp = dperp_fat*1e-9*ones([Ncomp 1]);
    w = ones([Ncomp 1]); w = w/sum(w);
    
elseif startsWith(system_name,'Watson_real_sticks_OP')
    OP = extractAfter(system_name,'Watson_real_sticks_OP_');
    load(fullfile(path,strcat('theta_phi_Watson_OP_', OP, '.mat')),'theta','phi');
    phi = phi'; % TO CORRECT
    mean_theta = theta;
    mean_phi = phi;
    [dpar, dperp, theta, phi, w] = set_of_real_sticks(mean_theta, mean_phi, path);
    
    
% Fibers and tumors systems

elseif startsWith(system_name,'composite_real_stick_and_ballf')
    f_ball = str2double(extractAfter(system_name,'composite_real_stick_and_ballf_'));
    [dpar, dperp, theta, phi, w] = get_real_stick(0, 0, path);
    D_iso = 3;
    dpar = [dpar ; D_iso*1e-9];
    dperp = [dperp ; D_iso*1e-9];
    theta = [theta ; 0];
    phi = [phi ; 0];
    w = [(1-f_ball)*w ; f_ball]; w = w/sum(w);
    
elseif startsWith(system_name,'composite_real_stick_and_CSFf')
    f_ball = str2double(extractAfter(system_name,'composite_real_stick_and_CSFf_'));
    [dpar, dperp, theta, phi, w] = get_real_stick(0, 0, path);
    D_iso_mean = 3*1e-9;
    sigma_value = 0.1*1e-9;
    D_iso = linspace(0.01*1e-9,3*1e-9,500)';
    theta_iso = zeros(size(D_iso));
    phi_iso = theta_iso;
    w_iso = exp(-1/2*(D_iso-D_iso_mean).^2/sigma_value^2); w_iso = w_iso/sum(w_iso);
    dpar = [dpar ; D_iso];
    dperp = [dperp ; D_iso];
    theta = [theta ; theta_iso];
    phi = [phi ; phi_iso];
    w = [(1-f_ball)*w ; f_ball*w_iso]; w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(D_iso, w_iso);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_edema.pdf'))
        clf(f)
    end
    
elseif startsWith(system_name,'composite_real_stick_and_edemaf')
    f_ball = str2double(extractAfter(system_name,'composite_real_stick_and_edemaf_'));
    [dpar, dperp, theta, phi, w] = get_real_stick(0, 0, path);
    D_iso_mean = 2.2*1e-9;
    sigma_value = 0.1*1e-9;
    D_iso = linspace(0.01*1e-9,3*1e-9,500)';
    theta_iso = zeros(size(D_iso));
    phi_iso = theta_iso;
    w_iso = exp(-1/2*(D_iso-D_iso_mean).^2/sigma_value^2); w_iso = w_iso/sum(w_iso);
    dpar = [dpar ; D_iso];
    dperp = [dperp ; D_iso];
    theta = [theta ; theta_iso];
    phi = [phi ; phi_iso];
    w = [(1-f_ball)*w ; f_ball*w_iso]; w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(D_iso, w_iso);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_edema.pdf'))
        clf(f)
    end
   
    
elseif startsWith(system_name,'composite_edema_and_gliomaf')
    f_tumor = str2double(extractAfter(system_name,'composite_edema_and_gliomaf_'));
    D_iso = linspace(0.01*1e-9,3*1e-9,500)';
    D_edema_mean = 2.2*1e-9;
    sigma_edema = 0.1*1e-9;
    D_glioma_mean = 0.7*1e-9;
    sigma_glioma = 0.1*1e-9;
    dpar = D_iso;
    dperp = D_iso;
    theta = zeros(size(D_iso));
    phi = theta;
    w = (1-f_tumor)*exp(-1/2*(D_iso-D_edema_mean).^2/sigma_edema^2) + f_tumor*exp(-1/2*(D_iso-D_glioma_mean).^2/sigma_glioma^2); 
    w = w/sum(w);
    %
    if is_plot
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        plot(D_iso, w);
        saveas(gcf,strcat(pwd, '/Distributions_systems/distribution_', system_name, '.pdf'))
        clf(f)
    end
    
  
    
end

end