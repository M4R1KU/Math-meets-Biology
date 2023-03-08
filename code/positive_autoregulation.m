function positive_autoregulation
    v_max = 4;
    K_diss = 10^(-16);
    k_dr = 0.2; 
    k_s = 5;
    k_dp = 0.7;

    tspan = [0 50];
    y0 = [0 0.01];
    [t,y] = ode45(@(t,y) odefcn(t, y, v_max, K_diss, k_dr, k_s, k_dp), tspan, y0);

    plot(t,y);
    legend('[RNA]','[P]');
end

function dydt = odefcn(t, y, v_max, K_diss, k_dr, k_s, k_dp)
  dydt = zeros(2,1);
  dydt(1) = v_max*(y(2)/(K_diss+y(2)))-k_dr*y(1);
  dydt(2) = k_s*y(1)-k_dp*y(2);
end