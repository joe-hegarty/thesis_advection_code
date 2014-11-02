function u_new = lax_wendroff(u, a, dx, dt)
    %LAX_WENDROFF Performs second-order lax_wendroff update

    % Assume u contains m+1 points from [0,1] inclusive
    m_plus_1 = length(u);
    m = m_plus_1 - 1;
    
    % Create array of m+1 points, the same size as u
    u_new = zeros(1, m_plus_1);
    
    % Assuming zero-based indexing for most of this algorithm, only process
    % indices 1,...,m.  Will adjust for one-based indexing at the end.
    for i = 1 : m;
        cn = a * dt / dx;
        
        i0 = i - 1;
        if i0 == 0;
            i0 = m;
        end
        
        i1 = i + 1;
        if i1 == m + 1;
            i1 = 1;
        end
        
        % Correct indexing to one-based indexing here only
        u_new(i+1) = u(i+1) - 0.5 * cn * (u(i1+1) - u(i0+1)) + 0.5 * cn * cn * (u(i1+1) - 2.0 * u(i+1) + u(i0+1));
    end
    
    % Set u(0) = u(1) using one-based indexing
    u_new(1) = u_new(m+1);
end

