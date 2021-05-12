close all
clear all

current_dir = '/Users/alexis_reymbaut/Desktop/Research_Lund/Codes/Matrix_Gamma/Test_dtd_mv_gamma';
Watson_path = [current_dir '/Watson'];

xps = load('xps_b0_6_10_21_46_extremes_mix_b0_lin_0_10_16_21_sph_2_2_2_0_pla_0_0_16_21.mat');
xps = xps.xps;

opt = get_opt;

SNR = 200;

% "unimodal_setvdison_0.8"
% "unimodal_setvdison_2.5"
% "unimodal_setvddeltan_0.75"
% "Watson_fat_structures_OP_0.01"
% "composite_real_stick_and_CSFf_0.5"

[dpar, dperp, theta, phi, w] = choose_simulated_system("composite_real_stick_and_CSFf_0.5", Watson_path);
dtd = dtd_par2dist(dpar,dperp,theta,phi,w);
m_true = dtd_dtd2m(dtd,opt);
dps_true = dtd_1d_fit2param(m_true);
s_true = dtd_1d_fit2data(m_true, xps);
s = sqrt((s_true + 1/SNR*randn([xps.n 1])).^2 + (1/SNR*randn([xps.n 1])).^2); % Adding Rician noise

method_name = 'dtd_mv_gamma';
tic
m = feval([method_name '_1d_data2fit'],s,xps,opt);
toc
dps = mio_1d_fit2param(method_name,m);
s_fit = feval([method_name '_1d_fit2data'],m,xps);
chisqn = msf_chisqn(s,s_fit);

method_name = 'dtd_covariance';
m_dtd = feval([method_name '_1d_data2fit'],s,xps,opt);
dps_dtd = mio_1d_fit2param(method_name,m_dtd);
s_fit_dtd = feval([method_name '_1d_fit2data'],m_dtd,xps);
chisqn = msf_chisqn(s,s_fit_dtd);

mdiso_true = dps_true.mdiso*1e9
mdiso_dtd_mv_gamma = dps.mdiso*1e9
% dps_dtd.mdiso*1e9
fprintf("\n")
vdiso_true = dps_true.vdiso*1e18
vdiso_dtd_mv_gamma = dps.vdiso*1e18
% dps_dtd.vdiso*1e18
fprintf("\n")
msdanison_true = dps_true.msdanison
msdanison_dtd_mv_gamma = dps.msdanison
% dps_dtd.msdanison

[~,ind_sorted] = sort(s,'descend');
x = 1:xps.n;


figure;
scatter(x,s(ind_sorted),20, 'filled', 'k')
hold on
scatter(x,s_fit(ind_sorted),20, 'b')
hold on
scatter(x,s_fit_dtd(ind_sorted),20, [0 0.8 0])