%% visualizing btotal(t)

clear
% define constants

r_m = 1560; % Europa radius (km)
B_prim = 250; % approximate background field amplitude nT
A = 0.95; % amplitude response factor
phi = 0.14; % phase lag in radians (around 8 deg)
synodic_period = 11.23 * 3600; % 11.23 hr period
omega = 2*pi / (synodic_period); % synodic freq. (rad/sec), 


% inducing field direction
e0x = [1;0;0];
e0y = [0;1;0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% create spatial grid 
grid_limit = 7 * r_m; % view up to 7 planetary radii out 
[X,Y,Z] = meshgrid(linspace(-grid_limit, grid_limit, 60), ...
    linspace(-grid_limit,grid_limit, 60), ...
    linspace(-grid_limit, grid_limit, 60));

% set up animation time array (60 steps over one full period)
num_frames = 60;
time_steps = linspace(0, synodic_period, num_frames);

% create figure once outside the loop
fig = figure('Color', 'b');
set(gcf, 'Position', [100, 100, 800, 600]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% animation loop 
for t = time_steps

    % compute secondary field components (bx,by,bz) at every point 
    % initialize arrays to hold magnetic field vector components
    Bx_spatial = zeros(size(X));
    By_spatial = zeros(size(X));
    Bz_spatial = zeros(size(X));
    
    % loop through each grid point to calculate the spatial part of eqn (3)
    for i = 1:numel(X)

        % define position vector r from center of moon
        r_vec = [X(i); Y(i); Z(i)];
        r = norm(r_vec);

        % mask field inside the moon to avoid division by 0
        if r < r_m
            continue;
        end


        % X EQUATION
        % calculate dot prod (r . e0x)
        r_dot_e0x = dot(r_vec, e0x);

        % evaluate spatial vector bracket: 3*(r . e0)*r - (r^2)*e0
        spatial_bracketx = 3 * r_dot_e0x * r_vec - (r^2) * e0x;

        % compute spatial coefficient: - (A * B_prim * r_m^3) / (2 * r^5)
        coefx = - (A * B_prim * r_m^3) * (cos (omega * t - phi)) / (2 *r^5);


        % Y EQUATION
        % calculate dot prod (r . e0y)
        r_dot_e0y = dot(r_vec, e0y);

        % evaluate spatial vector bracket: 3*(r . e0)*r - (r^2)*e0
        spatial_brackety = 3 * r_dot_e0y * r_vec - (r^2) * e0y;

        % compute spatial coefficient: - (A * B_prim * r_m^3) / (2 * r^5)
        coefy = - (A * B_prim * r_m^3) * (sin (omega * t - phi)) / (2 *r^5);


      



        % store vector components
        Bprim_vector = [B_prim * cos(omega * t); B_prim * sin(omega * t); 0];
        B_vector = (coefx * spatial_bracketx) + (coefy * spatial_brackety) + Bprim_vector;
        Bx_spatial(i) = B_vector(1);
        By_spatial(i) = B_vector(2);
        Bz_spatial(i) = B_vector(3);    

    end

    
    Bx = Bx_spatial;
    By = By_spatial;
    Bz = Bz_spatial;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    clf; % clear figure window so frames don't overlap



    hold on;


    % draw Europa as 3d sphere
    [sx, sy, sz] = sphere(24);
    mesh(sx*r_m, sy*r_m, sz*r_m, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor',[0.5 0.5 0.5]);
    alpha(0.4);

    % set seeds (dynamic)
    % find the current primary angle of the dipole field in x-y plane
    dipole_angle = atan2(sin(omega * t - phi), cos(omega * t - phi));

    % base ring template (orthogonal to dipole axis)
    theta = linspace(0, 2*pi, 16);
    seed_radius = 1.05 * r_m;

    % generate ring points centered around origin 
    % construct them relative to dipole axis
    local_long = seed_radius * 0.85;
    local_r = seed_radius * 0.5;

    % positive pole ring template
    x1_temp = ones(size(theta)) * local_long;
    y1_temp = local_r * cos(theta);
    z1 = local_r * sin(theta);

    % rotate template coordinated by dipole's current orientation angle
    seed_X1 = x1_temp * cos(dipole_angle) - y1_temp * sin(dipole_angle);
    seed_Y1 = x1_temp * sin(dipole_angle) + y1_temp * cos(dipole_angle);

    % negative pole ring template
    x2_temp = ones(size(theta)) * -local_long;
    y2_temp = local_r * cos(theta);
    z2 = local_r * sin(theta);

    seed_X2 = x2_temp * cos(dipole_angle) - y2_temp * sin(dipole_angle);
    seed_Y2 = x2_temp * sin(dipole_angle) + y2_temp * cos(dipole_angle);


    % combine everything into the master seed list
    seed_X = [seed_X1, seed_X2];
    seed_Y = [seed_Y1, seed_Y2];
    seed_Z = [z1, z2];

    % plot mag field vectors

    % calculate total magnetic field magnitude at every grid point
    B_mag = sqrt(Bx.^2 + By.^2 + Bz.^2);

    % get raw spatial vertices of the streamlines
    verts1 = stream3(X, Y, Z, Bx, By, Bz, seed_X(:), seed_Y(:), seed_Z(:));
    verts2 = stream3(X, Y, Z, -Bx, -By, -Bz, seed_X(:), seed_Y(:), seed_Z(:));
    all_verts = [verts1; verts2];

    % loop through each streamline and plot as colored patch
    for k = 1:numel(all_verts)
        v = all_verts{k};
        if isempty(v), continue; end

        % interpolate the field strength at each vertex along specific line
        line_colors = interp3(X, Y, Z, B_mag, v(:,1), v(:,2), v(:,3));

        % create multicolored line using patch
        patch('XData', v(:,1), 'YData', v(:,2), 'ZData', v(:,3), ...
            'FaceColor', 'none', 'EdgeColor', 'interp', ...
            'CData', line_colors, 'LineWidth', 1.5);

    end

    % apply colormap and colorbar scaling
    colormap(jet); % 'jet' or 'turbo' goes from blue (low) to red (high)
    c = colorbar;
    ylabel(c, 'Field Strength (nT)', 'FontSize', 11);

    % set color limits: because B drops off like 1/r^5, capping the max color 
    % at 1.3 * B_prim keeps the external fields from looking completely blue
    clim([0, B_prim * 1.3]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % grid formatting
    axis equal; 
    xlim([-grid_limit, grid_limit]);
    ylim([-grid_limit, grid_limit]);
    zlim([-grid_limit, grid_limit]);
    xlabel('x position (km)', 'FontSize', 11);
    ylabel('y position (km)', 'FontSize', 11);
    zlabel('z position (km)', 'FontSize', 11);

    % display current time in hours on the title
    title(sprintf('B-total @ t = %.2f hours', t/3600), 'FontSize', 12);
    grid on;
    box on;

    xlim([-4875 4832])
    ylim([-4835 4872])
    zlim([-4853 4853])


    view([-270 90])

    


    drawnow; % force matlab to draw graphics pipeline immediately
    pause(0.05); % execution speed

end