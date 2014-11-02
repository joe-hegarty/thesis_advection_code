function u_new = beam_warming(u, a, dx, dt)
    %BEAM_WARMING Performs second-order beam_warming update

    % Assume u contains m+1 points from [0,1] inclusive
    m_plus_1 = length(u);
    m = m_plus_1 - 1;
    
    % Create array of m+1 points, the same size as u
    u_new = zeros(1, m_plus_1);
    
    % Assuming zero-based indexing for most of this algorithm, only process
    % indices 1,...,m.  Will adjust for one-based indexing at the end.
    for i = 1 : m;
        cn = a * dt / dx;
        
        i_1 = i - 1;
        if i_1 == 0;
            i_1 = m;
        end
        
        i_2 = i_1 - 1;
        if i_2 == 0;
            i_2 = m;
        end
        
        % Correct indexing to one-based indexing here only
        u_new(i+1) = u(i+1) - 0.5 * cn * (3.0 * u(i+1) - 4.0 * u(i_1+1) + u(i_2+1)) + 0.5 * cn * cn * (u(i+1) - 2.0 * u(i_1+1) + u(i_2+1));
    end
    
    % Set u(0) = u(1) using one-based indexing
    u_new(1) = u_new(m+1);
end

