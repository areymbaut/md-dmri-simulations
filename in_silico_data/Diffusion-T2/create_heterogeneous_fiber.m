function structure_out = create_heterogeneous_fiber(structure_info)

%% Book-keeping
case_name = structure_info.case_name;
output_dir = structure_info.output_dir;
method = structure_info.method;

if ~isfolder(output_dir)
    msf_mkdir(output_dir);
end

%% Fiber properties
N = structure_info.N;
mean_diso = structure_info.mean_diso;
mean_ddelta = structure_info.mean_ddelta;
mean_theta = structure_info.mean_theta;
mean_phi = structure_info.mean_phi;
if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
    mean_t = structure_info.mean_t;
end
fraction = structure_info.fraction;
relative_std = structure_info.relative_std;
compartment_names = structure_info.compartment_names;
nb_compartments = length(mean_diso);

if nb_compartments == 1
    colors = [0 0 1];
elseif nb_compartments == 2
    colors = [[1 0 0]; [0 0 1]];
elseif nb_compartments == 3
    colors = [[1 0 0]; [0 0.8 0]; [0 0 1]];
end

%% Construct intra-fiber distributions (no dispersion for now)
diso = zeros(N, nb_compartments);
ddelta = zeros(N, nb_compartments);
if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
    t = zeros(N, nb_compartments);
end
theta = zeros(N, nb_compartments);
phi = zeros(N, nb_compartments);
w = zeros(N, nb_compartments);

for n = 1:nb_compartments
    diso_comp = mean_diso(n).*linspace((1-4*relative_std(n)), (1+4*relative_std(n)), N)';
    ddelta_comp = mean_ddelta(n).*linspace((1-4*relative_std(n)), (1+4*relative_std(n)), N)';
    if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
        t_comp = mean_t(n).*linspace((1-4*relative_std(n)), (1+4*relative_std(n)), N)';
    end
    theta_comp = mean_theta(n).*ones([N 1]);
    phi_comp = mean_phi(n).*ones([N 1]);
    
    dummy = linspace(-4*relative_std(n), 4*relative_std(n),N)';
    w_comp = exp(-dummy.^2./(2*relative_std(n)^2));
    w_comp = fraction(n).*w_comp./sum(w_comp);
    
    % Enforce the sign of covariances between quantities
    % C[diso, ddelta] < 0
    % C[diso, t] < 0
    % C[ddelta, t] > 0
    diso_comp = flip(diso_comp);
    
    diso(:,n) = diso_comp;
    ddelta(:,n) = ddelta_comp;
    theta(:,n) = theta_comp;
    phi(:,n) = phi_comp;
    if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
        t(:,n) = t_comp;
    end
    w(:,n) = w_comp;
end

% Test relative std
% n = 1;
% std_diso = sqrt(sum(diso(:,n).^2.*w(:,n)/fraction(n)) - sum(diso(:,n).*w(:,n)/fraction(n))^2)/sum(diso(:,n).*w(:,n)/fraction(n))
% std_ddelta = sqrt(sum(ddelta(:,n).^2.*w(:,n)/fraction(n)) - sum(ddelta(:,n).*w(:,n)/fraction(n))^2)/sum(ddelta(:,n).*w(:,n)/fraction(n))
% std_t = sqrt(sum(t(:,n).^2.*w(:,n)/fraction(n)) - sum(t(:,n).*w(:,n)/fraction(n))^2)/sum(t(:,n).*w(:,n)/fraction(n))

%% Create output structure
dpar = diso(:).*(1 + 2.*ddelta(:));
dperp = diso(:).*(1 - ddelta(:));
r = 1./t(:);
structure_out.dpar = dpar;
structure_out.dperp = dperp;
structure_out.theta = theta(:);
structure_out.phi = phi(:);
if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
    structure_out.r = r;
end
structure_out.w = w(:);

%% Plot
set(0, 'defaultLegendInterpreter','latex');
lw = 2;
global_font_size = 16;
global_font_size_labels = 17;

f = figure('Position', [0,0,1500,1500], 'Units', 'pixels', 'visible', 'off');
f.PaperOrientation = 'landscape';
f.PaperUnits = 'normalized';
f.PaperPosition = [0 0 1 1];
% set(gcf, 'Renderer', 'painters')

width = 0.4;
height = 0.2;
inter_v = 0.1;
offset_left = 0;
offset_bottom = 0.02;

x_left = offset_left + (1 - width)/2;

if strcmp(method, 'dtd')
    y_bottom = offset_bottom + (1 - 2*height - inter_v)/2;
    axh_diso = axes('position', [x_left y_bottom + height + inter_v width height]);
    axh_ddelta = axes('position', [x_left y_bottom width height]);
    hold(axh_diso, 'on')
    hold(axh_ddelta, 'on')
elseif strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
    y_bottom = offset_bottom + (1 - 3*height - 2*inter_v)/2;
    axh_diso = axes('position', [x_left y_bottom + 2*(height + inter_v) width height]);
    axh_ddelta = axes('position', [x_left y_bottom + (height + inter_v) width height]);
    axh_t = axes('position', [x_left y_bottom width height]);
    hold(axh_diso, 'on')
    hold(axh_ddelta, 'on')
    hold(axh_t, 'on')
end

for n = 1:nb_compartments
    diso_to_plot = squeeze(diso(:,n))*1e9;
    ddelta_to_plot = squeeze(ddelta(:,n));
    w_to_plot = squeeze(w(:,n));
    color = colors(n,:);
    
    [max_w, ind_max_w] = max(w_to_plot);
    
    plot(axh_diso, diso_to_plot, w_to_plot, '-', 'Linewidth', 2, 'color', color)
    plot(axh_ddelta, ddelta_to_plot, w_to_plot, '-', 'Linewidth', 2, 'color', color)
    
    if strcmp(method, 'dtr2d') 
        t_to_plot = squeeze(t(:,n))*1e3;
        plot(axh_t, t_to_plot, w_to_plot, '-', 'Linewidth', 2, 'color', color)
    elseif strcmp(method, 'dtr1d')
        t_to_plot = squeeze(t(:,n))*1e3;
        plot(axh_t, t_to_plot, w_to_plot, '-', 'Linewidth', 2, 'color', color)
    end
    
    text(axh_diso, diso_to_plot(ind_max_w), 1.01*max_w, compartment_names{n}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', global_font_size_labels, 'Color', color)
end

plot(axh_diso, (sum(diso(:).*w(:)))*1e9*[1 1], ylim(axh_diso), '--', 'Linewidth', 2, 'Color', 0.7*[1 1 1]);
plot(axh_ddelta, (sum(ddelta(:).*w(:)))*[1 1], ylim(axh_ddelta), '--', 'Linewidth', 2, 'Color', 0.7*[1 1 1]);
if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
    plot(axh_t, (sum(t(:).*w(:)))*1e3*[1 1], ylim(axh_t), '--', 'Linewidth', 2, 'Color', 0.7*[1 1 1]);
end

set(axh_diso, 'LineWidth', lw, 'TickDir', 'out', 'FontSize', global_font_size, 'TickLabelInterpreter', 'latex')
set(axh_ddelta, 'LineWidth', lw, 'TickDir', 'out', 'FontSize', global_font_size, 'TickLabelInterpreter', 'latex')
if strcmp(method, 'dtr2d') || strcmp(method, 'dtr1d')
    set(axh_t, 'LineWidth', lw, 'TickDir', 'out', 'FontSize', global_font_size, 'TickLabelInterpreter', 'latex')
end

xlabel(axh_diso, '$D_\mathrm{iso}$ ($\mu$m$^2$/ms)', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
xlabel(axh_ddelta, '$D_\Delta$', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
ylabel(axh_diso, '$P(D_\mathrm{iso})$', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
ylabel(axh_ddelta, '$P(D_\Delta)$', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
if strcmp(method, 'dtr2d')
    xlabel(axh_t, '$T_2$ (ms)', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
    ylabel(axh_t, '$P(T_2)$', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
elseif strcmp(method, 'dtr1d')
    xlabel(axh_t, '$T_1$ (s)', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
    ylabel(axh_t, '$P(T_1)$', 'Fontsize', global_font_size_labels, 'Interpreter', 'latex');
end
% set(xlab, 'Units', 'Normalized', 'Position', [0.5, -0.45, 0]);

xlim(axh_diso, [0.5, 1.5])
set(axh_diso, 'XTick', 0.5:0.25:1.5)
xlim(axh_ddelta, [0, 1])
if strcmp(method, 'dtr2d')
    xlim(axh_t, [40, 130])
elseif strcmp(method, 'dtr1d')
    xlim(axh_t, [0, 1])
end

saveas(f, fullfile(output_dir, strcat(case_name, '.pdf')))
clf(f)

end