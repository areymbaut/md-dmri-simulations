close all
clearvars

%% General parameters
method = 'dtd';
xps_name = 'xps_60points';
xps_file = fullfile(pwd, strcat(xps_name, '.mat'));

% SNRs and number of noise realizations for each SNR
SNR_list = [5 10 15 30 60 90 200];
N_SNR = 100; 

%% Load xps and create output directory
load(xps_file, 'xps');

output_dir = fullfile(pwd, xps_name);
if ~exist(output_dir, 'dir')
    msf_mkdir(output_dir);
end

% Suppress b0 signals
% ind = round(xps.b*1e-9,1) == 0;
% if nnz(ind) > 0
%     xps = mdm_xps_subsample(xps, ~ind);
%     struc.xps = xps;
%     save(xps_file, '-struct', 'struc');
%     load(xps_file, 'xps');
% end

%% Prepare system structure

% Optional system information to plot the ground-truth distributions
structure_info.do_plot = 1;
structure_info.case_name = '4_1iso_3way_crossing';
structure_info.compartment_names = {'iso'; 'aniso1'; 'aniso2'; 'aniso3'};
black = [0 0 0]; red = [1 0 0]; green = [0 0.8 0]; blue = [0 0 1];
structure_info.colors = [black; blue; green; red];

% Mandatory system information
N = 100;
structure_info.N = N;
structure_info.method = method;
structure_info.output_dir = output_dir;
structure_info.mean_diso = [2.5; 0.75; 0.8 ; 0.85]*1e-9;
structure_info.mean_ddelta = [0 0.9 0.85 0.8];
structure_info.mean_theta = [0 pi/2 pi/4 0];
structure_info.mean_phi = [0 0 pi/2 0];
structure_info.dispersion = [0 0 0 0];
structure_info.fraction = 0.25*ones([1, 4]);
structure_info.relative_std = 0.02.*ones(size(structure_info.mean_diso));

structure_out = create_heterogeneous_system(structure_info);
diso = structure_out.diso(:);
ddelta = structure_out.ddelta(:);
dpar = structure_out.dpar(:);
dperp = structure_out.dperp(:);
theta = structure_out.theta(:);
phi = structure_out.phi(:);
w = structure_out.w(:);
w = w/sum(w);

% Ground-truth statistical descriptors
sddelta = ddelta.^2;
mdiso = sum(w.*diso);
msddelta = sum(w.*sddelta);
vdiso = sum(w.*diso.^2) - mdiso^2;
vsddelta = sum(w.*sddelta.^2) - msddelta^2;
cvdisosddelta = sum(w.*(diso-mdiso).*(sddelta - msddelta));

%% Ground-truth signal
opt = mdm_opt();
opt = dtd_opt(opt);
opt.dtd.n_out = length(dpar); % Need to keep all the components
dtd = dtd_par2dist(dpar, dperp, theta, phi, w);
m = dtd_dtd2m(dtd, opt);
s_true = dtd_1d_fit2data(m, xps);

s_true = reshape(s_true, 1, 1, 1, xps.n);
mdm_nii_write(s_true, fullfile(output_dir, 'ground_truth_signal.nii.gz'));

%% Noise realizations
output_dir_signals = fullfile(output_dir, 'noisy_signals');
for it_SNR = 1:length(SNR_list)
    SNR = SNR_list(it_SNR);
    output_dir_signals_SNR = fullfile(output_dir_signals, strcat('SNR', num2str(SNR)));
   
    for n = 1:N_SNR
        s = sqrt((s_true + 1/SNR*randn(size(s_true))).^2 + (1/SNR*randn(size(s_true))).^2); % Adding Rician noise
        mdm_nii_write(s_true, fullfile(output_dir_signals_SNR, strcat('signal', num2str(n), '.nii.gz')));
    end
end