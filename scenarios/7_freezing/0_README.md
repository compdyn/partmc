
# Immersion freezing scenario

[reference]

This is a scenario demonstrating how to use PartMC to simulate the immersion freezing process of multiple species INPs.

There are a total of 4*2=8 simulation groups, consisting of four types of INP populations and two temperature curve scenarios. The four types of INPs are 100% illite, 100% Fe2O3, a 50% illite and 50% Fe2O3 external mixture, and an internal mixture. The two temperature curves are a constant -20 degrees Celsius and a steady cooling from -10 to -30 degrees Celsius. The simulation time for each is 10 minutes.

exp1: 100% illite, constant temperature
exp2: 100% Fe2O3, constant temperature
exp3: external mixture, constant temperature
exp3: internal mixture, constant temperature
exp5: 100% illite, steady cooling
exp6: 100% Fe2O3, steady cooling
exp7: external mixture, steady cooling
exp8: internal mixture, steady cooling

1. Run 1_run.sh. This is the shell scripts for running all simulations.
2. Run 2_process.sh. This is the shell scripts for extract data from netcdf files.
3. Run 3_draw.py. This creates the figure showing the frozen fraction time series in each simulation. (Reproduces the Figure 8 in Tang et al., 2025)
4. (Optional) Run 4_clean.sh. This deletes all files created by the processes above.
