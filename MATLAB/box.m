function y = box(x, s, e)
    %BOX Function going from 0 to 1 and then back to 0
    %   x coordinates at each index (x should be monotonic)
    %   s start of the box
    %   e end of the box
    %
    %   Assuming that the box function is periodic, i.e. repeats every unit
    y = zeros(1, length(x));
    for i = 2:length(x);
        if x(i) > s - floor(s) && (e - floor(e) < s - floor(s) || x(i) < e - floor(e));
            y(i) = 1.0;
        end
        if x(i) < e - floor(e) && (s - floor(s) > e - floor(e) || x(i) > s - floor(s));
            y(i) = 1.0;
        end
    end
    y(1) = y(length(y));
end
