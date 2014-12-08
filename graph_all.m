function graph_all(config_filename)
close all;

%Read Config File:
configID = fopen(config_filename);
delta_space = str2double(fgetl(configID));
grid_width = str2double(fgetl(configID));
grid_height = str2double(fgetl(configID));
birth_rate_filename = fgetl(configID);
t_filename = fgetl(configID);
i_filename = fgetl(configID);
v_filename = fgetl(configID);
tissue_type_filename = fgetl(configID);
delta_t = fgetl(configID);
number_of_timesteps = fgetl(configID);
result_t_filename = fgetl(configID);
result_i_filename = fgetl(configID);
result_v_filename = fgetl(configID);

fclose(configID);
%make all 6 plots (T,I,V, and the same afterwards - you know this from
%the config file reading)

%Now we read in all the data
T_ID = fopen(t_filename);
T_ID
grid_width
grid_height
plot_T = fread(T_ID, [grid_width, grid_height], 'double');
fclose(T_ID);

I_ID = fopen(i_filename);
plot_I = fread(I_ID, [grid_width, grid_height], 'double');
fclose(I_ID);

V_ID = fopen(v_filename);
plot_V = fread(V_ID, [grid_width, grid_height], 'double');
fclose(V_ID);

result_T_ID = fopen(result_t_filename);
result_plot_T = fread(result_T_ID, [grid_width, grid_height], 'double');
fclose(result_T_ID);

result_I_ID = fopen(result_i_filename);
result_plot_I = fread(result_I_ID, [grid_width, grid_height], 'double');
fclose(result_I_ID);

result_V_ID = fopen(result_v_filename);
result_plot_V = fread(result_V_ID, [grid_width, grid_height], 'double');
fclose(result_V_ID);

%And now we plot:
y = linspace(0,grid_height * delta_space,grid_height);
x = linspace(0,grid_width * delta_space,grid_width);

%Scaling down the size for very large grids just so matlab doesn't cry
if (grid_height * grid_width > 10000)
    plot_T = plot_T(1:10:end,1:10:end);
    plot_I = plot_I(1:10:end,1:10:end);
    plot_V = plot_V(1:10:end,1:10:end);
    result_plot_T = result_plot_T(1:10:end,1:10:end);
    result_plot_I = result_plot_I(1:10:end,1:10:end);
    result_plot_V = result_plot_V(1:10:end,1:10:end);
    
    x = x(1:10:end);
    y = y(1:10:end);
end
    
subplot(2,3,1);
surf(x,y,plot_T.');
set(gca,'xlim',[0 grid_width]);
set(gca,'ylim',[0 grid_height]);
title('starting T population');
shading interp
subplot(2,3,2);
surf(x,y,plot_I.');
set(gca,'xlim',[0 grid_width]);
set(gca,'ylim',[0 grid_height]);
title('starting I population');
shading interp
subplot(2,3,3);
surf(x,y,plot_V.');
set(gca,'xlim',[0 grid_width]);
set(gca,'ylim',[0 grid_height]);
title('starting V population');
shading interp
subplot(2,3,4);
surf(x,y,result_plot_T.');
set(gca,'xlim',[0 grid_width]);
set(gca,'ylim',[0 grid_height]);
title('ending T population');
shading interp
subplot(2,3,5);
surf(x,y,result_plot_I.');
set(gca,'xlim',[0 grid_width]);
set(gca,'ylim',[0 grid_height]);
title('ending I population');
shading interp
subplot(2,3,6);
surf(x,y,result_plot_V.');
set(gca,'xlim',[0 grid_width]);
set(gca,'ylim',[0 grid_height]);
title('ending V population');
shading interp