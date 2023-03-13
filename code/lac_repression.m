function lac_repression
    v_max = 0.4;
    K_diss = 10^(-6);
    k_dr = 0.36; 
    k_s = 0.4;
    k_dp = 0.056;

    K_diss_dna = 10^(-10);
    K_diss_lac_i = 6*10^(-6);


    k_lac = (3000/60)*10^(-4);

    k_dlactase = 0.056;

    lac_i = 10^(-8);

    tspan = [0 200];
    y0 = [0.005 0 0.001];
    [t,y] = ode45(@(t,y) odefcn(t, y, k_lac, v_max, K_diss_dna, K_diss_lac_i, k_dr, k_s, k_dlactase, lac_i), tspan, y0);

    plot(t,y);
    legend('[Lactose]','[Lactase]', '[mRNA]');
end

function dydt = odefcn(t, y, k_lac, v_max, K_diss_dna, K_diss_lac_i, k_dr, k_transl, k_dlactase, lac_i)
  dydt = zeros(3,1);
  dydt(1) = -k_lac*y(2);
  dydt(2) = k_transl*y(3)-k_dlactase*y(2);
  dydt(3) = v_max*(K_diss_dna/(K_diss_dna+(lac_i*(y(1)/(K_diss_lac_i+y(1))))))-k_dr*y(3);
end