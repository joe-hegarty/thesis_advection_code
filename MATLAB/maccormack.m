function u_new = maccormack(u, a, dx, dt)
    %MACCORMACK Performs a 2nd advection with limit clamping
    u_n_1 = semi_lagrangian(u, a, dx, dt);
    u_n = semi_lagrangian(u_n_1, a, dx, -dt);

    e = u_n - u;
    u_new = u_n_1 - 0.5 * e;
    
    %u_new = limiter_revert(u, a, dx, dt, u_new, u_n_1);    
    %u_new = limiter_clamp(u, a, dx, dt, u_new);
end
