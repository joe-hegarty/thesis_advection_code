function y = step( x, a )
    %STEP Function going from 1 to 0 a units from the origin
    %   x coordinates at each index (x should be monotonic)
    %   a specifies signed distance from the origin of step
    y = zeros(1, length(x));
    for i = 1:length(x);
        if x(i) < a;
            y(i) = 1.0;
        else
            y(i) = 0.0;
        end
    end
end
