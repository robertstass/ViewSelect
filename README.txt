Select cryo-EM particles from a viewing direction distribution plot.

This script will take an exported cryosparc .csg or .cs file and display a viewing direction distribution plot.
The user can then draw circles on the plot to select particles that are of a particular orientation.
Separate .csg/.cs files containing these particles are then saved that can be imported back into cryosparc.

Note that if you supply a .csg file, the corresponding .cs file must be present in the same directory.
New .csg files are only output if a .csg file is input originally (Recommended as this makes re-importing easier).

Usage:

Start by exporting particles from the output tab of a 3D refinement job (Actions -> Export). Then:

python view_select.py cryosparc_exported_file.csg
or select the file from the GUI simply with
python view_select.py

For help:
python view_select.py -h

Dependencies:
cryosparc-tools (install with "pip install cryosparc-tools")
matplotlib (Tested with v3.7.0. Some earlier versions don't have the ability to rotate an ellipse selection).

Installation:
Go to the "View Select" directory and use the following commands:
conda env create -f view_select.yml
conda activate view_select
pip install -r requirements.txt
(optionally) export PATH=<path_to>/ViewSelect:$PATH

Then it can be run using
"conda active view_select" followed by "python view_select.py"

Why do this?
For small and/or featureless samples, it is common to get one or two views that dominate all others resulting
in a preferred orientation problem. It is then optimal to try to enrich other views but these can sometimes be
hard to find at the 2D classification stage. The 3D refinement stage can sometimes help but it can be difficult
to relate hotspots in a viewing direction distribution plot back to the 2D classes. There is a partial solution
to this which is to generate a .bild file for viewing in chimera using pyem (https://github.com/asarnow/pyem).
However, "View Select" allows the user to select particles from the plot that can then be used for further 2D classification
and/or particle picking with tools like Topaz.

Written by Robert Stass, Bowden group, STRUBI/OPIC (2024)
