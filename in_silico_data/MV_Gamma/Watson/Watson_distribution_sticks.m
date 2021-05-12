function [theta,phi] = Watson_distribution_sticks(kappa, do_plot)

close all

if nargin < 2
    do_plot = 0;
end

%% Generate the Watson-distributed points
N_sticks = 100;
points = randWatson(N_sticks, [0,0,1], kappa);
x = points(:,1);
y = points(:,2);
z = points(:,3);
[theta,phi] = cartesian2spherical_unit_sphere(x,y,z);

%% Distribution visualization
if do_plot
    ax = polaraxes;
    polarscatter(ax, phi, sin(theta))
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
end
