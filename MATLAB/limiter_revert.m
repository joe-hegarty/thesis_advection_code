function u_new = limiter_revert(u, a, dx, dt, u_1, u_1_hat)
    %LIMITER_REVERT If the new value is greater than or less than the values
    %               of the two interpolated sample points, revert to standard
    %               semi-langrangian
    
    % Assume u contains m+1 points from [0,1] inclusive
    m_plus_1 = length(u);
    m = m_plus_1 - 1;
    
    u_new = u_1;
    
    % Assuming zero-based indexing for most of this algorithm, only process
    % indices 1,...,m.  Will adjust for one-based indexing at the end.
    for i = 1 : m;
        x = i - a * dt / dx;
        
        while x < 0;
            x = x + m;
        end
        while x > m;
            x = x - m;
        end
        
        i0 = floor(x);
        i1 = i0 + 1;
        if i0 == 0;
            i0 = m;
        end
        if i1 == m + 1;
            i1 = 1;
        end
        
        % Correct indexing to one-based indexing here only
        u_min = min(u(i0+1), u(i1+1));
        u_max = max(u(i0+1), u(i1+1));
        if u_new(i+1) > u_max || u_new(i+1) < u_min;
            u_new(i+1) = u_1_hat(i+1);
        end
    end
    
    % Set u(0) = u(1) using one-based indexing
    u_new(1) = u_new(m+1);
end
