clear all

in_vivo_data_path = fullfile(pwd,'in_vivo_data');
data_name = 'BRAIN_FWF_MERGED';
xps_name = 'BRAIN_FWF_MERGED_xps';
output_path = fullfile(in_vivo_data_path, 'in_vivo_results_dtd_mv_gamma_do_multiple_s0');
nb_slice = 14;

%% Connect to data
[S, nii_h] = mdm_nii_read(fullfile(in_vivo_data_path,strcat(data_name,'.nii.gz')));
S = S(:,:,nb_slice,:);
mdm_nii_write(S, fullfile(in_vivo_data_path,strcat(data_name, '_slice',num2str(nb_slice),'.nii.gz')), nii_h);
s.nii_fn = fullfile(in_vivo_data_path,strcat(data_name, '_slice',num2str(nb_slice),'.nii.gz'));
% s.mask_fn = fullfile(i, 'data_mask.nii.gz');

%% xps and opt
s.xps = mdm_xps_load(fullfile(in_vivo_data_path, xps_name));
opt = mdm_opt();
opt = dtd_mv_gamma_opt(opt);
opt.do_mask = 1;
opt.do_data2fit = 1;
opt.do_fit2param = 1;
opt.do_param2nii = 1;
opt.do_multiple_s0 = 1; % CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

msf_mkdir(output_path);
paths.pdf_path = fullfile(output_path, 'pdf_maps');
paths.nii_path = fullfile(output_path, 'nii_maps');
paths.mfs_fn = fullfile(output_path, 'mfs.mat');
paths.chisq_fn = fullfile(output_path, 'chisq.mat');
paths.ind_fn = fullfile(output_path, 'ind.mat');
paths.dps_fn = fullfile(output_path, 'dps.mat');
paths = mdm_paths(paths);

tic
if opt.do_data2fit && isempty(gcp('nocreate'))
    parpool(4)
end
% parpool(3, 'IdleTimeout', Inf)
nii_fn = dtd_mv_gamma_pipe(s, paths, opt);
toc