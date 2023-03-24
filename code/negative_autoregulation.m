function negative_autoregulation
    v_max = 15; % nano M min^-1
    K_diss = 10^(-2); % nano M

    k_s = 20; % min^-1
    
    k_dr = 0.23; % min^-1 (Quelle: Bistable Behaviour of the lac operon)
    k_dp = 0.011; % min^-1

    tspan = [0 100];
    y0 = [0 0];
    [t,y] = ode45(@(t,y) negative_diff_system(y, v_max, K_diss, k_dr, k_s, k_dp), tspan, y0);

    plot(t,y(:,1));
    ylabel('[RNA] / nM');

    yyaxis right
    plot(t,y(:,2));
    ylabel('[P] / nM');

    xlabel('Zeit / min')
    legend({'[RNA]', '[P]'}, 'Location', 'east');
    legend('boxoff');
    axis padded;
end

function dydt = negative_diff_system(y, v_max, K_diss, k_dr, k_s, k_dp)
  dydt = zeros(2,1);
  dydt(1) = v_max*K_diss/(K_diss+y(2))-k_dr*y(1);
  dydt(2) = k_s*y(1)-k_dp*y(2);
end