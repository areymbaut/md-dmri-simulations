function opt = get_opt

opt = mdm_opt();
opt = dtd_gamma_opt(opt);
opt = dtd_mv_gamma_opt(opt);
opt = dtd_pa_opt(opt);
opt = dtd_covariance_opt(opt);
opt = dtd_opt(opt);