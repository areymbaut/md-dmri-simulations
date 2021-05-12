clear
close all

set(0, 'defaultLegendInterpreter','latex');

model = 'dtd';
opt = mdm_opt();

% Prepare paths
data_path = fullfile(pwd,'in_vivo_data');
mfs_path = fullfile(data_path,'in_vivo_results_dtd_mv_gamma_mc_do_multiple_s0');
% paths.voxel_path     = fullfile(data_path, analysis, 'voxel');

% Connect to data
nii_path = fullfile(mfs_path,'nii_maps');
s.nii_fn = fullfile(nii_path, 'BRAIN_FWF_MERGED_mc_slice14.nii.gz');
s.mask_fn = fullfile(nii_path, 'BRAIN_FWF_MERGED_mc_slice14_mask.nii.gz');
s.xps = mdm_xps_load(fullfile(nii_path, 'BRAIN_FWF_MERGED_mc_xps.mat'));
xps = s.xps;

%% Sort Signal
b_r = sqrt(xps.u(:,1).^2 + xps.u(:,2).^2 + xps.u(:,3).^2);
b_theta = acos(xps.u(:,3)./b_r);
b_phi = atan2(xps.u(:,2),xps.u(:,1));

acq_mtrx = [round(xps.b/1e9,4) round(xps.b_delta,4) b_theta b_phi];
[acq_mtrx_sort,indx_sort] = sortrows(acq_mtrx);


%% Voxels of interest
voxel_WM = [46 39 1];
voxel_GM = [61 19 1];
voxel_CSF = [44 46 1];
voxel_crossing = [36 62 1];

mfs = mdm_mfs_load(fullfile(mfs_path, 'mfs.mat'));
measured_signal = mdm_nii_read(s.nii_fn);

measured_signal_WM = squeeze(measured_signal(voxel_WM(1),voxel_WM(2),voxel_WM(3),:));
measured_signal_GM = squeeze(measured_signal(voxel_GM(1),voxel_GM(2),voxel_GM(3),:));
measured_signal_CSF = squeeze(measured_signal(voxel_CSF(1),voxel_CSF(2),voxel_CSF(3),:));
measured_signal_crossing = squeeze(measured_signal(voxel_crossing(1),voxel_crossing(2),voxel_crossing(3),:));
fitted_signal_WM = dtd_mv_gamma_1d_fit2data(squeeze(mfs.m(voxel_WM(1),voxel_WM(2),voxel_WM(3),:)), xps);
fitted_signal_GM = dtd_mv_gamma_1d_fit2data(squeeze(mfs.m(voxel_GM(1),voxel_GM(2),voxel_GM(3),:)), xps);
fitted_signal_CSF = dtd_mv_gamma_1d_fit2data(squeeze(mfs.m(voxel_CSF(1),voxel_CSF(2),voxel_CSF(3),:)), xps);
fitted_signal_crossing = dtd_mv_gamma_1d_fit2data(squeeze(mfs.m(voxel_crossing(1),voxel_crossing(2),voxel_crossing(3),:)), xps);

measured_signal_WM = double(measured_signal_WM(indx_sort));
measured_signal_GM = double(measured_signal_GM(indx_sort));
measured_signal_CSF = double(measured_signal_CSF(indx_sort));
measured_signal_crossing = double(measured_signal_crossing(indx_sort));
fitted_signal_WM = double(fitted_signal_WM(indx_sort));
fitted_signal_GM = double(fitted_signal_GM(indx_sort));
fitted_signal_CSF = double(fitted_signal_CSF(indx_sort));
fitted_signal_crossing = double(fitted_signal_crossing(indx_sort));

%% Plot Sorted acquisition scheme

% Reset coordinates and sizes
lw_contour = .5;
lw_axes = 1;
lw_plot = .5;
fs_axes = 5.5;
papersize = [8.7 8.7];
sz_marker = .5;

width = .6;
inter_v = 0.04;
height = .07;
tick_l = .03*[1 1];

shift_value = .08;

left = (1-width)/2;
bottom = (1-8*height-7*inter_v)/2+0.03;
y8 = bottom;
y7 = y8 + height + inter_v;
y6 = y7 + height + inter_v;
y5 = y6 + height + inter_v;
y4 = y5 + height + inter_v;
y3 = y4 + height + inter_v;
y2 = y3 + height + inter_v;
y1 = y2 + height + inter_v;

figure('visible', 'off');
clf

axh_2 = axes('position',[left y1 width height]);
%     set(axh,'ColorOrder',[1 0 0; 0 .7 0; 0 0 1],'NextPlot', 'replacechildren');
ph_2 = plot(1:xps.n,xps.b(indx_sort)/1e9,'.');
%     ph = scatter(1:numel(signal),log10(signal),scale*mean(a),'filled');
shift = shift_value*(2-0);
axis([-.08 1.05*xps.n -shift 2+shift])
set(ph_2,'MarkerSize',sz_marker,'MarkerEdgeColor',[0 0 0])%,'MarkerFaceColor',[0 0 0])
set(axh_2,'LineWidth',lw_axes,'FontSize',fs_axes)
set(axh_2,'XTick',[0 100 200 300],'XTickLabel',[],'YTick',[0 1 2],'TickDir','out','TickLength',tick_l,'Box','off')
axh_2.YAxis.TickLabelInterpreter = 'latex';
ylab = ylabel(axh_2, '$b\;(\mathrm{ms}/\mu\mathrm{m}^2)$', 'FontSize',fs_axes, 'Interpreter', 'latex');
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);

axh_3 = axes('position',[left y2 width height]);
%     set(axh,'ColorOrder',[1 0 0; 0 .7 0; 0 0 1],'NextPlot', 'replacechildren');
ph_3 = plot(1:xps.n,xps.b_delta(indx_sort),'.');
%     ph = scatter(1:numel(signal),log10(signal),scale*mean(a),'filled');
shift = shift_value*(1-(-.5));
axis([-.08 1.05*xps.n -0.5-shift 1+shift])
set(ph_3,'MarkerSize',sz_marker,'MarkerEdgeColor',[0 0 0])%,'MarkerFaceColor',[0 0 0])
set(axh_3,'LineWidth',lw_axes,'FontSize',fs_axes)
set(axh_3,'XTick',[0 100 200 300],'XTickLabel',[],'YTick',[-.5 0 1],'TickDir','out','TickLength',tick_l,'Box','off')
axh_3.YAxis.TickLabelInterpreter = 'latex';
ylab = ylabel(axh_3, '$b_\Delta$', 'FontSize',fs_axes, 'Interpreter', 'latex');
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);

axh_4 = axes('position',[left y3 width height]);
%     set(axh,'ColorOrder',[1 0 0; 0 .7 0; 0 0 1],'NextPlot', 'replacechildren');
ph_4 = plot(1:xps.n,b_theta(indx_sort),'.');
%     ph = scatter(1:numel(signal),log10(signal),scale*mean(a),'filled');
shift = shift_value*(pi-(0));
axis([-.08 1.05*xps.n 0-shift pi+shift])
set(ph_4,'MarkerSize',sz_marker,'MarkerEdgeColor',[0 0 0])%,'MarkerFaceColor',[0 0 0])
set(axh_4,'LineWidth',lw_axes,'FontSize',fs_axes)
axh_4.YAxis.TickLabelInterpreter = 'latex';
ylab = ylabel(axh_4, '$\Theta$', 'FontSize',fs_axes, 'Interpreter', 'latex');
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(axh_4,'XTick',[0 100 200 300],'XTickLabel',[],'YTick',[0 pi/2 pi],'YTickLabel',{'0', '$\pi/2$', '$\pi$'},'TickDir','out','TickLength',tick_l,'Box','off')

axh_5 = axes('position',[left y4 width height]);
%     set(axh,'ColorOrder',[1 0 0; 0 .7 0; 0 0 1],'NextPlot', 'replacechildren');
ph_5 = plot(1:xps.n,abs(b_phi(indx_sort)-pi),'.');
%     ph = scatter(1:numel(signal),log10(signal),scale*mean(a),'filled');
shift = shift_value*(2*pi-(0));
axis([-.08 1.05*xps.n 0-shift 2*pi+shift])
set(ph_5,'MarkerSize',sz_marker,'MarkerEdgeColor',[0 0 0])%,'MarkerFaceColor',[0 0 0])
set(axh_5,'LineWidth',lw_axes,'FontSize',fs_axes)
axh_5.XAxis.TickLabelInterpreter = 'latex';
axh_5.YAxis.TickLabelInterpreter = 'latex';
ylab = ylabel(axh_5, '$\Phi$', 'FontSize',fs_axes, 'Interpreter', 'latex');
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(axh_5,'XTick',[0 100 200 300],'XTickLabel',[],'YTick',[0 pi 2*pi],'YTickLabel',{'0', '$\pi$', '$2\pi$'},'TickDir','out','TickLength',tick_l,'Box','off')

axh_6 = axes('position',[left y5 width height]);
hold(axh_6,'on')
%     set(axh,'ColorOrder',[1 0 0; 0 .7 0; 0 0 1],'NextPlot', 'replacechildren');
ph_6 = plot(axh_6, 1:xps.n, measured_signal_CSF/max(measured_signal_CSF), '.');
ph_fit_6 = plot(axh_6, 1:xps.n, fitted_signal_CSF/max(fitted_signal_CSF), '.');
set(ph_6, 'MarkerSize', sz_marker, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])%,'MarkerFaceColor',[0 0 0])
set(ph_fit_6, 'MarkerSize', 0.5*sz_marker, 'MarkerFaceColor', [floor(50/255) floor(125/255) 1], 'MarkerEdgeColor', [floor(50/255) floor(125/255) 1])%,'MarkerFaceColor',[0 0 0])
shift = shift_value*(1-(0));
axis([-.08 1.05*xps.n 0-shift 1+shift])
set(axh_6,'LineWidth',lw_axes,'FontSize',fs_axes)
axh_6.XAxis.TickLabelInterpreter = 'latex';
axh_6.YAxis.TickLabelInterpreter = 'latex';
ylab = ylabel(axh_6, '$\tilde{S}_\mathrm{CSF}$', 'FontSize',fs_axes, 'Interpreter', 'latex');
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(axh_6,'XTick',[0 100 200 300],'XTickLabel',[],'YTick',[0 0.5 1],'YTickLabel',{'0', '0.5', '1'},'TickDir','out','TickLength',tick_l,'Box','off')

axh_7 = axes('position',[left y6 width height]);
hold(axh_7,'on')
%     set(axh,'ColorOrder',[1 0 0; 0 .7 0; 0 0 1],'NextPlot', 'replacechildren');
ph_7 = plot(axh_7, 1:xps.n, measured_signal_GM/max(measured_signal_GM), '.');
ph_fit_7 = plot(axh_7, 1:xps.n, fitted_signal_GM/max(fitted_signal_GM), '.');
set(ph_7, 'MarkerSize', sz_marker, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])%,'MarkerFaceColor',[0 0 0])
set(ph_fit_7, 'MarkerSize', 0.5*sz_marker, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0])%,'MarkerFaceColor',[0 0 0])
shift = shift_value*(1-(0));
axis([-.08 1.05*xps.n 0-shift 1+shift])
set(axh_7,'LineWidth',lw_axes,'FontSize',fs_axes)
axh_7.XAxis.TickLabelInterpreter = 'latex';
axh_7.YAxis.TickLabelInterpreter = 'latex';
ylab = ylabel(axh_7, '$\tilde{S}_\mathrm{GM}$', 'FontSize',fs_axes, 'Interpreter', 'latex');
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(axh_7,'XTick',[0 100 200 300],'XTickLabel',[],'YTick',[0 0.5 1],'YTickLabel',{'0', '0.5', '1'},'TickDir','out','TickLength',tick_l,'Box','off')

axh_1 = axes('position',[left y7 width height]);
hold(axh_1,'on')
%     set(axh,'ColorOrder',[1 0 0; 0 .7 0; 0 0 1],'NextPlot', 'replacechildren');
ph_1 = plot(axh_1, 1:xps.n, measured_signal_WM/max(measured_signal_WM), '.');
ph_fit_1 = plot(axh_1, 1:xps.n, fitted_signal_WM/max(fitted_signal_WM), '.');
set(ph_1, 'MarkerSize', sz_marker, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])%,'MarkerFaceColor',[0 0 0])
set(ph_fit_1, 'MarkerSize', 0.5*sz_marker, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0])%,'MarkerFaceColor',[0 0 0])
shift = shift_value*(1-(0));
axis([-.08 1.05*xps.n 0-shift 1+shift])
set(axh_1,'LineWidth',lw_axes,'FontSize',fs_axes)
axh_1.XAxis.TickLabelInterpreter = 'latex';
axh_1.YAxis.TickLabelInterpreter = 'latex';
ylab = ylabel(axh_1, '$\tilde{S}_\mathrm{WM}$', 'FontSize',fs_axes, 'Interpreter', 'latex');
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(axh_1,'XTick',[0 100 200 300],'XTickLabel',[],'YTick',[0 0.5 1],'YTickLabel',{'0', '0.5', '1'},'TickDir','out','TickLength',tick_l,'Box','off')

axh_8 = axes('position',[left y8 width height]);
hold(axh_8,'on')
%     set(axh,'ColorOrder',[1 0 0; 0 .7 0; 0 0 1],'NextPlot', 'replacechildren');
ph_8 = plot(axh_8, 1:xps.n, measured_signal_crossing/max(measured_signal_crossing), '.');
ph_fit_8 = plot(axh_8, 1:xps.n, fitted_signal_crossing/max(fitted_signal_crossing), '.');
set(ph_8, 'MarkerSize', sz_marker, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])%,'MarkerFaceColor',[0 0 0])
set(ph_fit_8, 'MarkerSize', 0.5*sz_marker, 'MarkerFaceColor', [1 0.65 0], 'MarkerEdgeColor', [1 0.65 0])%,'MarkerFaceColor',[0 0 0])
shift = shift_value*(1-(0));
axis([-.08 1.05*xps.n 0-shift 1+shift])
set(axh_8,'LineWidth',lw_axes,'FontSize',fs_axes)
axh_8.XAxis.TickLabelInterpreter = 'latex';
axh_8.YAxis.TickLabelInterpreter = 'latex';
ylab = ylabel(axh_8, '$\tilde{S}_\mathrm{crossing}$', 'FontSize',fs_axes, 'Interpreter', 'latex');
xlab = xlabel(axh_8, '$n_\mathrm{acq}$', 'FontSize',fs_axes, 'Interpreter', 'latex');
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(xlab, 'Units', 'Normalized', 'Position', [0.5, -0.8, 0]);
set(axh_8,'XTick',[0 100 200 300],'YTick',[0 0.5 1],'YTickLabel',{'0', '0.5', '1'},'TickDir','out','TickLength',tick_l,'Box','off')

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize);
saveas(gcf, fullfile(pwd,'acq_ordered_ISMRM2020.pdf'))


% %% Plot 5D Sampling scheme
% 
% % bpar = xps.b.*(1 + 2*xps.b_delta);
% % bperp = xps.b.*(1 - xps.b_delta);
% 
% x = xps.b/1e9; x = x(:); x(~isfinite(x)) = 0;
% y = xps.b_delta; y = y(:); y(~isfinite(y)) = 0;
% z = xps.tr;  z = z(:); z(~isfinite(z)) = 0;
% a = 15;
% 
% acpro.x = x;
% acpro.y = y;
% acpro.z = z;
% acpro.a = a*ones(xps.n,1);
% 
% acpro.r = abs(xps.u(:,1)); acpro.r = acpro.r(:);
% acpro.g = abs(xps.u(:,2)); acpro.g = acpro.g(:);
% acpro.b =abs(xps.u(:,3)); acpro.b = acpro.b(:);
% acpro.bright = 0*ones(xps.n,1);
% 
% 
% %%%%%%% Plot Acquisition Scheme as a Scatter plot
% 
% % Coordinates and sizes
% lw_contour = 1;
% lw_axes = 2.5;
% lw_plot = 1;
% fs_axes = 10*7/4.2;
% tick_l = 0.05*[1 1];
% papersize = [10 10];
% left = .15;
% bottom = .15;
% width = .8;
% height = .8;
% 
% figure('visible', 'off'); 
% clf
% 
% axh = axes('position',[left bottom width height]);
% %     ph = scatter3(x,y,z,a,[0 0 0]);
% contourpars.Nx = 75; contourpars.Ny = contourpars.Nx; contourpars.Nz = contourpars.Nx; contourpars.Nlevels = 20;
% axpars.xmin = -0.25; axpars.xmax = 2.25; axpars.ymin = -.75; axpars.ymax = 1.25; axpars.zmin = 0; axpars.zmax = 5.6;
% 
% [axh,ph,hcontour] = dist_3d_scattercontourplot(axh,acpro,contourpars,axpars);
% %     view(30,30)
% %     xmin = 0; xmax = 2; ymin = -.5; ymax = 1; zmin = 0; zmax = .151;
% %     axis([xmin xmax ymin ymax zmin zmax])
% set(ph,'LineWidth',lw_plot,'MarkerFaceColor',[0 0 0])
% set(hcontour,'LineWidth',.5*lw_plot)%,'Color',.75*[1 0 0])
% set(axh,'LineWidth',lw_axes,'FontSize',fs_axes)
% axis(axh,'square')
% set(axh,'XTick',0:.5:2,'YTick',-.5:.5:1,'ZTick',0:1:5,'TickLength',tick_l)
% set(axh,'XTickLabel',{'0','','1','','2'},'YTickLabel',{'-0.5','0','0.5','1'},'ZTickLabel',{'0','1','2','3','4','5'})
% set(axh,'Projection','perspective')
% axh.XAxis.TickLabelInterpreter = 'latex';
% axh.YAxis.TickLabelInterpreter = 'latex';
% axh.ZAxis.TickLabelInterpreter = 'latex';
% 
% fig = gcf;
% fig.Color = 'none';
% fig.Renderer = 'Painters';
% fig.InvertHardcopy = 'off';
% set(fig, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize);
% saveas(gcf, fullfile(pwd,'3d_acq_scheme.pdf'))


