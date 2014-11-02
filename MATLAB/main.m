% Start from a clean slate
clear; close all; clc;


% Experiment specific setup
title_name = 'BFECC';
method_name = 'bfecc';
init_func = 'box';
video_file = strcat('C:/Thesis_Media/Thesis_Videos/', method_name);
image_file = strcat('C:/Thesis_Media/Thesis_Images/', method_name);
output_video = false;
output_image = true;
close_figure = true;


% Figure configuration
label_font_size = 14;


% Video configuration
fps = 30;
slow_down = 4.0;
frame_mod = int32(1);


% Experiment configuration
time_total = 4.0;
m = 256;
a = 1.0;
cfl = 3.8;
over_sample = 16;

sin_scaling = 4;
a_0 = 0.225;
a_1 = 0.275;
b_0 = 0.7;
b_1 = 0.8;


% Derived values
dx = 1 / m;
x = 0 : dx : 1;
dt = abs(cfl * dx / a);
n_steps = int32(time_total / dt);
x_os = 0 : (dx / over_sample) : 1;
u = zeros(1, length(x));
u_os = zeros(1, length(x_os));

% Adjust video settings to get correct playback
if ((1 / dt) / slow_down) < fps;
    fps = (1 / dt) / slow_down;
else
    frame_mod = int32((1 / dt) / (slow_down * fps));
end


% Initialise u
if strcmp(init_func, 'box');
    u = box(x, a_0, a_1);
    u = u + box(x, b_0, b_1);

    u_os = box(x_os, a_0, a_1);
    u_os = u_os + box(x_os, b_0, b_1);
elseif strcmp(init_func, 'sin');
    u = sin_wave(x, 0, sin_scaling*2*pi, 0.5, 0.5);
    u_os = sin_wave(x_os, 0, sin_scaling*2*pi, 0.5, 0.5);
elseif strcmp(init_func, 'nyquist');
end


U = u;


% Initialise figure
screen_size = get(0, 'ScreenSize');
hFigure = figure('OuterPosition', [10 10 ((3*screen_size(3))/4) ((3*screen_size(4))/4)]);
hAxes = axes;
hDataExactU = plot(hAxes, x_os, u_os, 'Color', [0.8 0.8 0.8], 'LineWidth', 2);
axis(hAxes,[0 1 -0.2 1.2]);
title(title_name, 'FontSize', label_font_size);
xlabel('x', 'FontSize', label_font_size);
ylabel('u', 'FontSize', label_font_size);
hold on;

hDataU = plot(x, U, 'b', 'LineWidth', 2);

% Video output setup
vidObj = VideoWriter(video_file, 'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = fps;
if output_video;
    open(vidObj); %#ok<*UNRCH>
end

for i = 0:n_steps;
    % Update original values
    if strcmp(init_func, 'box');
        u_os = box(x_os, a_0 + double(i) * dt * a, a_1 + double(i) * dt * a);
        u_os = u_os + box(x_os, b_0 + double(i) * dt * a, b_1 + double(i) * dt * a);    
    elseif strcmp(init_func, 'sin');
        u_os = sin_wave(x_os, -double(i) * dt * a * sin_scaling*2*pi, sin_scaling*2*pi, 0.5, 0.5);
    end

    % Update results
    set(hDataExactU, 'XData', x_os, ...
                     'YData', u_os);
               
    set(hDataU, 'XData', x, ...
                'YData', U);

    if output_video;
        if mod(i, frame_mod) == 0 || i == n_steps;
            writeVideo(vidObj, getframe(hFigure));
        end
    end
    
    drawnow;

    % Advect
    if strcmp(method_name, 'forward_euler');
        U = forward_euler(U, a, dx, dt);
    elseif strcmp(method_name, 'lax_friedrichs');
        U = lax_friedrichs(U, a, dx, dt);
    elseif strcmp(method_name, 'upwind_euler');
        U = upwind_euler(U, a, dx, dt);
    elseif strcmp(method_name, 'semi_lagrangian');
        U = semi_lagrangian(U, a, dx, dt);
    elseif strcmp(method_name, 'lax_wendroff');
        U = lax_wendroff(U, a, dx, dt);
    elseif strcmp(method_name, 'beam_warming');
        U = beam_warming(U, a, dx, dt);
    elseif strcmp(method_name, 'maccormack');
        U = maccormack(U, a, dx, dt);
    elseif strcmp(method_name, 'bfecc');
        U = bfecc(U, a, dx, dt);
    elseif strcmp(method_name, 'minmod_slope');
        U = minmod_slope(U, a, dx, dt);
    elseif strcmp(method_name, 'superbee_slope');
        U = superbee_slope(U, a, dx, dt);
    end

end

if output_image;
    saveas(hFigure, image_file, 'png');
end

if output_video;
    close(vidObj);
end

if close_figure;
    close(hFigure);
end

fprintf('Experiment finished successfully!\n');
