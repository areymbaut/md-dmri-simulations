function [dpar, dperp, theta, phi, w] = get_real_stick_dispersive(mean_theta, mean_phi, OP_char, path)
    mean_dpar = 1.77*1e-9;
    mean_dperp = 0.31*1e-9;
    sigma = 0.09*1e-9;
    width_dpar = 0.4*1e-9;
    width_dperp = 0.14*1e-9;
    
    load(fullfile(path,strcat('theta_phi_Watson_OP_', OP_char, '.mat')),'theta','phi');
    phi = phi'; % TO CORRECT
    for i = 1:length(theta)
        z = cos(theta(i));
        if z < 0
            theta(i) = theta(i) - pi;
            phi(i) = phi(i) + pi;
        end
    end

    % Sort radially
    [~,ind_sort_angle] = sort(abs(sin(theta)));
    theta = theta(ind_sort_angle);
    phi = phi(ind_sort_angle);
    N = length(theta);
    
    if mean_theta ~= 0 || mean_phi ~= 0
        Ry = [[cos(mean_theta), 0, sin(mean_theta)];
            [0 1 0];
            [-sin(mean_theta), 0, cos(mean_theta)]];
        
        Rz = [[cos(mean_phi), -sin(mean_phi), 0];
            [sin(mean_phi), cos(mean_phi), 0];
            [0, 0, 1]];
        
        for i = 1:N
            cart_v = [sin(theta(i))*cos(phi(i)); sin(theta(i))*sin(phi(i)); cos(theta(i))];
            rot_cart_v = Rz*Ry*cart_v;
            [theta(i),phi(i)] = cartesian2spherical_unit_sphere(rot_cart_v(1), rot_cart_v(2), rot_cart_v(3));
        end
    end
    
    dpar = linspace(mean_dpar-width_dpar,mean_dpar+width_dpar,N)';
    dperp = linspace(mean_dperp-width_dperp,mean_dperp+width_dperp,N)';   
    [~,ind_sort_dpar] = sort(abs(dpar-mean_dpar));
    [~,ind_sort_dperp] = sort(abs(dperp-mean_dperp));
    dpar = dpar(ind_sort_dpar);
    dperp = dperp(ind_sort_dperp);
    
    w = exp(-1/2*(dpar-mean_dpar).^2/sigma^2); w = w/sum(w);
    w = sort(w,'descend');
   
%     diso = (dpar + 2*dperp)/3;
%     var(diso)/mean(diso)^2
    
%     mean_dpar = sum(w.*dpar);
%     std_dpar = sqrt(sum(w.*(dpar-mean_dpar).^2));
%     mean_dperp = sum(w.*dperp);
%     std_dperp = sqrt(sum(w.*(dperp-mean_dperp).^2));
%     
%     mean_dpar_value = mean_dpar*1e9
%     relative_error_dpar = std_dpar/mean_dpar
%     mean_dperp_value = mean_dperp*1e9
%     relative_error_dperp = std_dperp/mean_dperp
     
%     ind_w5percent = w>0.005;
%     cos_beta = cos(theta(ind_w5percent)).*cos(mean_theta) + sin(theta(ind_w5percent)).*sin(mean_theta).*cos(phi(ind_w5percent)-mean_phi);
%     mean_angular_dispersion = rad2deg(mean(abs(sum(w(ind_w5percent).*acos(cos_beta))/sum(w(ind_w5percent)))))
%     OP_w5percent = (3*sum(w(ind_w5percent).*cos_beta.^2)/sum(w(ind_w5percent))-1)/2
%     
%     cla
%     ax = polaraxes;
%     for i = 1:length(phi)
%         polarscatter(ax, phi(i), sin(theta(i)), 'b', 'filled', 'MarkerFaceAlpha', w(i)/max(w))
%         hold on
%         polarscatter(ax, phi(i), sin(theta(i)), 'r', 'filled', 'MarkerFaceAlpha', (1-w(i)/max(w))/10)
%         hold on
%     end
%     ax.ThetaDir = 'clockwise';
%     ax.ThetaZeroLocation = 'top';
    
    
    % '0.1', 0.6736 -
    % '0.3', 0.7981 -
    % '0.5', 0.9018 -
    % '0.8', 0.9565 -
    % '0.9', 0.9837 -
    
    
    % '0.2', 0.7647
    % '0.4', 0.8802
    % '0.6', 0.9362
    % '0.7', 0.9427
    
    