
# Immersion freezing scenario

Reference:
Tang, W. (2024). Particle-resolved simulations of immersion freezing with multi-species ice nucleating particles. (Master’s thesis, University of Illinois at Urbana-Champaign). https://hdl.handle.net/2142/124611

This is a scenario demonstrating how to use PartMC to simulate the immersion freezing process of multiple species INPs.

Species: Fe2O3 and illite (50:50)
Mixing state: internal mixture
Temperature: -10 to -30 degrees Celsius.
Immersion freezing scheme: ABIFM.

1. Run 1_run.sh. This is the shell scripts for running simulation.
2. Run 2_process.sh. This is the shell scripts for extract data from netcdf files.
3. Run 3_draw.py. This creates the figure showing the frozen fraction time series in each simulation. 
   (the output pdf file is out/TSs.pdf)
4. (Optional) Run 4_clean.sh. This deletes all files created by the processes above.
