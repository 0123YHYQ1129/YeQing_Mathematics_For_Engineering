function main
    % Define the ODE function
    odefun = @(t, y) [0, 1; -(1 + 2*cos(t)), 0] * y;
    
    % Set the time span and initial conditions
    tspan = [0, 2*pi];
    y0_1 = [1; 0];
    y0_2 = [0; 1];
    
    % Set options for ode45 with specified tolerances
    options = odeset('AbsTol', 1e-7, 'RelTol', 1e-8);

    % Call ode45 to solve the ODE with specified options
    [T1, Y1] = ode45(odefun, tspan, y0_1, options);
    [T2, Y2] = ode45(odefun, tspan, y0_2, options);

    % Construct the monodromy matrix
    monodromy_matrix = [Y1(end, :)', Y2(end, :)'];

    % Calculate eigenvalues
    eigenvalues = eig(monodromy_matrix);
    
    % Calculate stability
    magnitudes = abs(eigenvalues);
    disp('Magnitudes of the eigenvalues:');
    disp(magnitudes)
    if any(magnitudes > 1)
        disp('The system is unstable');
    else
        disp('The system is stable');
    end
    
    % Calculate λ
    lambda = (1 / (2 * pi)) * log(eigenvalues);
    disp('Eigenvalues λ:');
    disp(lambda);
end
