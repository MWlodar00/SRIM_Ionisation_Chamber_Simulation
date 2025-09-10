The following packages and programs need to be installed to run the script
- SRIM (installation folder needs to be specified in the config file)
- pysrim
- nyclyr
- periodictable
- pandas

The current (10/09/2025) version of pysrim does not allow for proper compound correction. To circumvent, the densities of different layers are scaled by a value obtained from material tables in the standalone SRIM program.
This scales the strengths of both the nuclear and electronic forces (only the latter should be concerned). However, since the reaction cross-section at the energy range (several MeV/u) is dominated by Coulomb scattering,
the difference between the approaches is minimal. The screenshot below shows the comparison between the results obtained with the standalone version and this code using a beam of Sn108 at 8MeV/u. In the future, the pysrim
package may be updated to include this feature.

<img width="1814" height="1064" alt="image" src="https://github.com/user-attachments/assets/798eee28-4c1d-4730-979f-5d272d42da64" />
