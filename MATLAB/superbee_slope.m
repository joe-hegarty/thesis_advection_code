function u_new = superbee_slope(u, a, dx, dt)
    %SUPERBEE_SLOPE Performs slope limited second-order update

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
        
        % superbee_slope
        sig_i = 0;
        sig1_i = 0;
        left_arg = (2.0*(u(i + 1) - u(i_1 + 1))) / dx;
        right_arg = (u(i1 + 1) - u(i + 1)) / dx;
        if left_arg * right_arg > 0;
            if abs(left_arg) < abs(right_arg);
                sig1_i = left_arg;
            else
                sig1_i = right_arg;
            end
        end
        sig2_i = 0;
        left_arg = (u(i + 1) - u(i_1 + 1)) / dx;
        right_arg = (2.0*(u(i1 + 1) - u(i + 1))) / dx;
        if left_arg * right_arg > 0;
            if abs(left_arg) < abs(right_arg);
                sig2_i = left_arg;
            else
                sig2_i = right_arg;
            end
        end
        if sig1_i * sig2_i > 0;
            if abs(sig1_i) < abs(sig2_i);
                sig_i = sig2_i;
            else
                sig_i = sig1_i;
            end
        end
        
        sig_i_1 = 0;
        sig1_i_1 = 0;
        left_arg = (2.0*(u(i_1 + 1) - u(i_2 + 1))) / dx;
        right_arg = (u(i + 1) - u(i_1 + 1)) / dx;
        if left_arg * right_arg > 0;
            if abs(left_arg) < abs(right_arg);
                sig1_i_1 = left_arg;
            else
                sig1_i_1 = right_arg;
            end
        end
        sig2_i_1 = 0;
        left_arg = (u(i_1 + 1) - u(i_2 + 1)) / dx;
        right_arg = (2.0*(u(i + 1) - u(i_1 + 1))) / dx;
        if left_arg * right_arg > 0;
            if abs(left_arg) < abs(right_arg);
                sig2_i_1 = left_arg;
            else
                sig2_i_1 = right_arg;
            end
        end
        if sig1_i_1 * sig2_i_1 > 0;
            if abs(sig1_i_1) < abs(sig2_i_1);
                sig_i_1 = sig2_i_1;
            else
                sig_i_1 = sig1_i_1;
            end
        end

        % Correct indexing to one-based indexing here only
        u_new(i + 1) = u(i + 1) - cn * (u(i + 1) - u(i_1 + 1)) - 0.5 * cn * (dx - a * dt) * (sig_i - sig_i_1);
    end
    
    % Set u(0) = u(1) using one-based indexing
    u_new(1) = u_new(m+1);
end

