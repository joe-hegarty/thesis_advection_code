function u_new = minmod_slope(u, a, dx, dt)
    %MINMOD_SLOPE Performs slope limited second-order update

    % Assume u contains m+1 points from [0,1] inclusive
    m_plus_1 = length(u);
    m = m_plus_1 - 1;
    
    % Create array of m+1 points, the same size as u
    u_new = zeros(1, m_plus_1);
    
    cn = a * dt / dx;

    % Assuming zero-based indexing for most of this algorithm, only process
    % indices 1,...,m.  Will adjust for one-based indexing at the end.
    for i = 1 : m;    
        i_1 = i - 1;
        if i_1 == 0;
            i_1 = m;
        end
        
        i_2 = i_1 - 1;
        if i_2 == 0;
            i_2 = m;
        end
        
        i1 = i + 1;
        if i1 == m + 1;
            i1 = 1;
        end
        
%         % Upwind Euler
%         sig_i = 0;
%         sig_i_1 = 0;
%         
%         % Lax-Wendroff
%         sig_i = (u(i1 + 1) - u(i + 1)) / dx;
%         sig_i_1 = (u(i + 1) - u(i_1 + 1)) / dx;

        % minmod_slope
        sig_i = 0;
        left_arg = (u(i + 1) - u(i_1 + 1)) / dx;
        right_arg = (u(i1 + 1) - u(i + 1)) / dx;
        if left_arg * right_arg > 0;
            if abs(left_arg) < abs(right_arg);
                sig_i = left_arg;
            else
                sig_i = right_arg;
            end
        end
        
        sig_i_1 = 0;
        left_arg = (u(i_1 + 1) - u(i_2 + 1)) / dx;
        right_arg = (u(i + 1) - u(i_1 + 1)) / dx;
        if left_arg * right_arg > 0;
            if abs(left_arg) < abs(right_arg);
                sig_i_1 = left_arg;
            else
                sig_i_1 = right_arg;
            end
        end

        % Correct indexing to one-based indexing here only
        u_new(i + 1) = u(i + 1) - cn * (u(i + 1) - u(i_1 + 1)) - 0.5 * cn * (dx - a * dt) * (sig_i - sig_i_1);
    end
    
    % Set u(0) = u(1) using one-based indexing
    u_new(1) = u_new(m+1);
end

