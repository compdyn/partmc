
# Urban plume scenario

N. Riemer, M. West, R. A. Zaveri, and R. C. Easter (2009) Simulating the evolution of soot mixing state with a particle-resolved aerosol model, _J. Geophys. Res._ 114(D09202), DOI: <http://dx.doi.org/10.1029/2008JD011073>.

This scenario simulates an air parcel moving over a polluted urban area during the daytime, and then moving up into an isolated boundary layer during the night. It has a realistic set of emissions and includes full chemistry and coagulation.

The processing scripts demonstrate two different ways of extracting data from PartMC: (1) using the `extract_*` programs to extract basic text-format data from the output NetCDF files, and (2) using a custom `urban_plume_process.F90` program to extract more complex data such as mixing state entropies.
