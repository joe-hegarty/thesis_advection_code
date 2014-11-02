function y = sin_wave(x, s, b, a, r)
    %SIN_WAVE Function going from 1 to 0 a units from the origin
    %   x coordinates at each index (x should be monotonic)
    %   s shift sin
    %   b bias the wavelength
    %   a amplitude
    %   r rise above the x-axis
    y = zeros(1, length(x));
    for i = 1:length(x);
        y(i) = a * sin(b * x(i) + s) + r;
    end
end
