function Watson_distribution_sticks_OP

kappa_test = 7;
OP_aim = 0.75;
count_while = 0;
while count_while < 1000
    [theta,phi] = Watson_distribution_sticks(kappa_test);
    OP = (3*mean(cos(theta).^2)-1)/2;
    if round(OP,3) == OP_aim
        mat_file.theta = theta;
        mat_file.phi = phi;
        save(['theta_phi_Watson_OP_' num2str(round(OP,1)) '.mat'],'-struct','mat_file');
        OP
        break
    end
    count_while = count_while + 1;
    round(OP,3)
end