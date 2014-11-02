function u_new = bfecc(u, v, dx, dt)
    %BFECC Performs a 2nd advection with limit clamping
    u_n_1 = semi_lagrangian(u, v, dx, dt);
    u_n = semi_lagrangian(u_n_1, v, dx, -dt);

    e = u_n - u;
    u_n = u - 0.5 * e;
    u_new = semi_lagrangian(u_n, v, dx, dt);
    
    %u_new = limiter_revert(u, v, dx, dt, u_new, u_n_1);    
    %u_new = limiter_clamp(u, v, dx, dt, u_new);
end
