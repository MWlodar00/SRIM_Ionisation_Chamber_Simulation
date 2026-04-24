This script only works on Windows.

The following packages and programs need to be installed to run the script
- SRIM (installation folder needs to be specified in the config file)
- pysrim
- nyclyr
- periodictable
- pandas

The current (10/09/2025) version of pysrim does not allow for proper compound correction. To circumvent, the densities of different layers are scaled by a value obtained from material tables in the standalone SRIM program. This scales the strengths of both the nuclear and electronic forces (only the latter should be concerned). However, since the reaction cross-section at the energy range (several MeV/u) is dominated by electrostatic forces, the difference between the approaches is minimal.

The code serves two main purposes:
1. Simulation of the ion range to gauge the ideal experimental conditions. When the *error_analysis* variable is set to *True* in the config file, the code includes the recommended 4% error on the stopping powers (https://www.srim.org/SRIM/SRIM2011.htm). Below is a comparison between simulation and experimental data for a Sn108 beam (with In contamination) at 8MeV/u and 92mbar of gas obtained during IS686. The agreement is more than satisfactory, but in most cases, one can expect a larger deviation.
<img width="1907" height="670" alt="image" src="https://github.com/user-attachments/assets/3badbaf6-2fe4-4265-8ad6-1468a5c448dc" />
2. A rough check of whether a particular beam mixture is resolvable. Using the past data on the resolution (obtained with an imperfect setup, so not of optimal quality), the code looks for a segment with the highest chances of resolving the beam components and plots what the final spectrum may look like. The resolution is quite hard to predict since simulating the energy straggling is very challenging, and the segment-to-segment variability in terms of electronics is non-negligible. Here's a comparison between simulated and recorded spectra for Sn108 (again). The resolution estimation worked quite well, but the spacing between the peaks is narrower than measured by about 1 MeV, hence the difference.
<img width="1853" height="568" alt="image" src="https://github.com/user-attachments/assets/8d52c856-898d-4210-af9c-57c9a7bab2d7" />
