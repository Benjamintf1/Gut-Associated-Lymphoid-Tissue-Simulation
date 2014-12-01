Gut-Associated-Lymphoid-Tissue-Simulation
=========================================

This is a simulation of Hiv or other viruses in gut associated lymphoid tissue with a focus on paralizing for large data sets.

To Compile: //TODO

To Run: //TODO

Configuration Details: //TODO: complete this section

A configuration file should be a file with one double, integer, or string value per line in the following order:
* (double) delta_space
* (integer) grid_width
* (integer) grid_height
* (string) birth_rate_filename:
* (string) t_filename:
* (string) i_filename:
* (string) v_filename:
* (string) tissue_type_filename
* (double) delta_t:
* (integer) number_of_timesteps:
* (string) result_t_filename:
* (string) result_i_filename:
* (string) result_v_filename:
* (integer) number_of_tissue_types
* (string) sub_config_for_tissue_type_0
* ...
* (string) sub_config_for_tissue_type_n

And for each tissue type, a sub-configuration file (with a name matching what is listed in the main configuration file) should be a file with one double or integer per line in the following order (specialized for the tissue type):
* (double) diffusion_tcells:
* (double) diffusion_infected:
* (double) diffusion_virus:
* (double) death_tcells:
* (double) death_infected:
* (double) death_virus:
* (integer) burst_rate:
* (double) transmission_vt:
* (double) transmission_it:


