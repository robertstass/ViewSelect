Select cryo-EM particles from a viewing direction distribution plot.

This script will take an exported cryosparc .csg/.cs file or a relion star file and display a viewing direction
distribution plot. The user can then draw circles on the plot to select particles that are of a particular orientation.
Separate .csg/.cs/.star files containing these particles are then saved that can be imported back into cryosparc/relion.

Note that if you supply a .csg file, the corresponding .cs file must be present in the same directory.
New .csg files are only output if a .csg file is input originally (Recommended as this makes re-importing easier).

![Alt text](example_plot.png?raw=true "Example plot")

Usage:\
For cryosparc, start by exporting particles from the output tab of a 3D refinement job (Actions -> Export).\
For relion, simply locate an appropriate star file from a 3D refinement.\
Then:\
python view_select.py exported_file.csg\
or select the file from the GUI simply with\
python view_select.py

For help:\
python view_select.py -h

Dependencies:\
cryosparc-tools (for cryosparc users. install with "pip install cryosparc-tools")\
starfile (by Alister Burt. For relion users. install with "pip install starfile")\
matplotlib (Tested with v3.7.0. Some earlier versions don't have the ability to rotate an ellipse selection).

Installation:\
Go to the "View Select" directory and use the following commands:\
conda env create -f view_select.yml\
conda activate view_select\
pip install -r requirements.txt\
(optionally) export PATH=<path_to>/ViewSelect:$PATH

Then it can be run using\
"conda activate view_select" followed by "python view_select.py"

On some HPC systems, the script can be sluggish due to high metadata traffic. It may therefore be better to containerize
it using something like apptainer/singularity/docker. Below is an example using the included .def file.\
apptainer build view_select.sif view_select_apptainer.def\
Then add the following to your ~/.bashrc file.\
alias view_select.py="apptainer exec -B /exafs1,/well /path/to/ViewSelect/view_select.sif python /path/to/ViewSelect/view_select.py"\
Then restart the terminal or do source ~/.bashrc\
Then to run simply use:\
view_select.py\
(Remember to change /path/to to the relevant full path. The -B /exafs1,/well argument is to mount additional drives to
the container. These will be specific to your system so edit as appropriate if needed.)

Why do this?\
For small and/or featureless samples, it is common to get one or two views that dominate all others resulting
in a preferred orientation problem. It is then optimal to try to enrich other views but these can sometimes be
hard to find at the 2D classification stage. The 3D refinement stage can sometimes help but it can be difficult
to relate hotspots in a viewing direction distribution plot back to the 2D classes. There is a partial solution
to this which is to generate a .bild file for viewing in chimera using pyem (https://github.com/asarnow/pyem).
Also native .bild file creation has now been added to cryosparc v4.4+. However, "View Select" allows the user to
select particles from the plot that can then be used for further 2D classification and/or particle picking with tools
like Topaz.

Further details:\
To use, draw one or more circles on the plot. If you click "New group" you can then draw more circles that are ultimately
saved into separate files. Use this to select different views. If two groups overlap each other, the overlapping particles
default to lower numbered groups. A separate file is saved which contains the remaining unselected particles. Click finish
to save the output.

If you click finish without any selections, an image of the plot is saved. This is an easy way to make cryosparc style
viewing direction distribution plots from relion star files.

Re-import .csg files into cryosparc with the Import Result Group job.

The values of rot, tilt and psi refer to euler angles in the relion convention, which are sequential rotations around the
z, y and z axis respectively. These euler angles have the effect of rotating the 3D map before generating the viewing
direction plot. It is helpful for selecting hotspots that have become distorted at the top and bottom of the plot.

The colors used can be customized for accessibility or aesthetic reasons using the --plot_cmap and --group_cmap command
line options. --cmap_limit can be set to apply a constant scale to the plot colormaps.

Written by Robert Stass, Bowden group, STRUBI/OPIC (2024)
