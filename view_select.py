#!/usr/bin/env python
print('Importing modules...')
import os
import sys
import math
import numpy as np
from threading import Thread, Event
import queue
import argparse
import copy
import time

from matplotlib import use as Use
Use("TkAgg")

from matplotlib import MatplotlibDeprecationWarning
from matplotlib import pyplot as plt
from matplotlib.widgets import EllipseSelector, RadioButtons
from matplotlib.transforms import Transform
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
from matplotlib.lines import Line2D
try:
    from matplotlib import colormaps as cm
except AttributeError:
    from matplotlib import cm as cm


import tkinter as tk
from tkinter import scrolledtext
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog

import warnings
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)

try:
    #from cryosparc.tools import CryoSPARC
    from cryosparc.dataset import Dataset
    import yaml
    cryosparc_support = True
except ImportError as e:
    cryosparc_support = False
    print(e)
    print('cryosparc-tools and PyYAML modules required for cryosparc support. Use: pip install cryosparc-tools')

try:
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    import starfile
    from starfile.writer import StarWriter
    relion_support = True
except ImportError as e:
    relion_support = False
    print(e)
    print('starfile python module required for relion star file support. Use: pip install starfile')

if not relion_support and not cryosparc_support:
    error_msg = 'Must have at least cryosparc-tools/pyYAML or starfile modules installed for cryosparc or relion support respectively.'
    raise ImportError(error_msg)







class ArgumentParser():

    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''
            Select cryo-EM particles from a viewing direction distribution plot.

            This script will take an exported cryosparc .csg/.cs file or a relion star file and display a viewing direction distribution
            plot. The user can then draw circles on the plot to select particles that are of a particular orientation.
            Separate .csg/.cs/.star files containing these particles are then saved that can be imported back into cryosparc/relion.
            
            Note that if you supply a .csg file, the corresponding .cs file must be present in the same directory.
            New .csg files are only output if a .csg file is input originally (Recommended as this makes re-importing easier).
            
            Usage:
            
            For cryosparc, start by exporting particles from the output tab of a 3D refinement job (Actions -> Export).
            For relion, simply locate an appropriate star file from a 3D refinement.
            Then:
            
            python view_select.py exported_file.csg
            or select the file from the GUI simply with
            python view_select.py
            
            For help:
            python view_select.py -h
            
            Dependencies:
            cryosparc-tools (for cryosparc users. install with "pip install cryosparc-tools")
            starfile (by Alister Burt. For relion users. install with "pip install starfile")
            matplotlib (Tested with v3.7.0. Some earlier versions don't have the ability to rotate an ellipse selection).
            
            Installation:
            Go to the "View Select" directory and use the following commands:
            conda env create -f view_select.yml
            conda activate view_select
            pip install -r requirements.txt
            (optionally) export PATH=<path_to>/ViewSelect:$PATH
            
            Then it can be run using
            "conda activate view_select" followed by "python view_select.py"
            
            Why do this?
            For small and/or featureless samples, it is common to get one or two views that dominate all others resulting
            in a preferred orientation problem. It is then optimal to try to enrich other views but these can sometimes be
            hard to find at the 2D classification stage. The 3D refinement stage can sometimes help but it can be difficult
            to relate hotspots in a viewing direction distribution plot back to the 2D classes. There is a partial solution
            to this which is to generate a .bild file for viewing in chimera using pyem (https://github.com/asarnow/pyem).
            However, "View Select" allows the user to select particles from the plot that can then be used for further 2D classification
            and/or particle picking with tools like Topaz.
            
            Further details:
            To use, draw one or more circles on the plot. If you click "New group" you can then draw more circles that are ultimately
            saved into separate files. Use this to select different views. If two groups overlap each other, the overlapping particles
            default to lower numbered groups. A separate file is saved which contains the remaining unselected particles. Click finish
            to save the output.
            
            If you click finish without any selections, an image of the plot is saved. This is an easy way to make cryosparc style viewing
            direction distribution plots from relion star files.
            
            Re-import .csg files into cryosparc with the Import Result Group job.
            
            The values of rot, tilt and psi refer to euler angles in the relion convention, which are sequential rotations around the
            z, y and z axis respectively. These euler angles have the effect of rotating the 3D map before generating the viewing
            direction plot. It is helpful for selecting hotspots that have become distorted at the top and bottom of the plot.
            
            The colors used can be customized for accessibility or aesthetic reasons using the --plot_cmap and --group_cmap command
            line options.
            
            Written by Robert Stass, Bowden group, STRUBI/OPIC (2024)
            ''')
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        add('--dont_save_image', action='store_true', help='Stop the script outputting an image at the end.')
        add('--plot_cmap', default=default_plot_cmap, help='Matplotlib colormap for the hexbin plot. (default: %(default)s)')
        add('--group_cmap', default=default_group_cmap, help='Matplotlib colormap for the selection group colors. (default: %(default)s)')
        add('input_dataset_path', nargs='?', default=None, help='File path of an exported cryosparc .csg/.cs file or a relion star file. Can leave blank to open a file in the GUI.')

        '''
        if len(sys.argv)==1:
                self.usage()
                sys.exit()
        '''
    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("Error: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

    def validate(self, args):
        if not args.input_dataset_path == None:
            if not os.path.exists(args.input_dataset_path):
                self.error("Input file '%s' not found." % args.input_dataset_path)
        try:
            cm.get_cmap(args.plot_cmap)
        except ValueError:
            self.error("Cannot load colormap with name '%s'." % args.plot_cmap)
        try:
            cm.get_cmap(args.group_cmap)
        except ValueError:
            self.error("Cannot load colormap with name '%s'." % args.group_cmap)



#Defaults
default_plot_cmap = 'jet'
default_group_cmap = 'tab10'


######### Main script ##########


class CryoDatafile:
    def __init__(self, dataset_path=None, cryosparc_support=True, relion_support=True):
        self.group_path = None #.csg
        self.dataset_path = dataset_path #.cs
        self.num_particles = None
        self.cryosparc_support = cryosparc_support
        self.relion_support = relion_support
        self.is_relion = False
        self.is_cryosparc = False
        self.rotation_matrices = None
        self.additional_dataset_files = []
        self.additional_dataset_output_paths = []

    def load_dataset(self):
        splitpath = os.path.splitext(self.dataset_path)
        if splitpath[1] == '.star':
            if not self.relion_support:
                error_msg = 'Install starfile module to load star files.'
                print_text(error_msg, color='red')
                raise ImportError(error_msg)
            self.is_relion = True
            self.is_cryosparc = False
            particles = self.load_rln_dataset()
        elif splitpath[1] == '.csg' or splitpath[1] == '.cs':
            if not self.cryosparc_support:
                error_msg = 'Install cryosparc-tools/pyYAML modules to load .csg/.cs files.'
                print_text(error_msg, color='red')
                raise ImportError(error_msg)
            self.is_relion = False
            self.is_cryosparc = True
            particles = self.load_cs_dataset()
        else:
            error_msg = 'File with extension "%s" not recognised.' % splitpath[1]
            print_text(error_msg, color='red')
            raise FileNotFoundError(error_msg)
        return particles

    def load_cs_dataset(self):
        splitpath = os.path.splitext(self.dataset_path)
        if splitpath[1] == '.csg':
            print_text('.csg file supplied.')
            print_text('Reading .csg file...')
            self.group_path = self.dataset_path
            self.load_csg()
        else:
            print_text('.cs file supplied. (no .csg file will be produced).')
        print_text('Reading .cs file...')
        self.particles = self.load_cs(self.dataset_path)
        self.num_particles = len(self.particles)
        return self.particles

    def load_rln_dataset(self):
        print_text('Reading .star file...')
        df = starfile.read(self.dataset_path)
        if 'optics' in df.keys():
            self.optics = df['optics']
        else:
            self.optics = None
        self.particles = df['particles']
        self.num_particles = len(self.particles)
        return self.particles

    def load_csg(self):
        with open(self.group_path, 'r') as file:
            self.group_data = yaml.safe_load(file)
        cs = self.group_data.get('results', {}).get('alignments3D', {}).get('metafile')
        num_items = self.group_data.get('results', {}).get('alignments3D', {}).get('num_items')
        if cs == None:
            error_msg = 'Error: .csg file must contain an "alignments3D" field.'
            print_text(error_msg, color='red')
            raise OSError(error_msg)
        if cs is not None and cs.startswith('>'):
            cs = cs[1:]
        dirname = os.path.dirname(self.group_path)
        if dirname != '':
            cs = os.path.join(dirname, cs)
        self.additional_dataset_files = []
        results_section = self.group_data.get('results', {})
        for result_label, result_value in results_section.items():
            metafile = result_value['metafile']
            metafile = metafile[1:] if metafile is not None and metafile.startswith('>') else metafile
            if metafile != cs:
                if result_value['num_items'] == num_items: #must be same as original otherwise ignore.
                    self.additional_dataset_files.append(metafile)
        self.additional_dataset_files = list(set(self.additional_dataset_files))
        for file_path in [cs] + self.additional_dataset_files:
            if not os.path.exists(file_path):
                error_msg = 'Error: When supplying a .csg file, the corresponding .cs file (%s) must be present in the same directory.' % (cs)
                print_text(error_msg, color='red')
                raise OSError(error_msg)
        self.dataset_path = cs

    def load_cs(self, file_path):
        return Dataset.load(file_path)

    def update_csg(self, new_metafile, new_num_items):
        #make a copy
        group_data = copy.deepcopy(self.group_data)
        old_metafile = os.path.basename(self.dataset_path)
        results_section = group_data.get('results', {})
        for result_label, result_value in results_section.items():
            if 'metafile' in result_value:
                if result_value['metafile'] == '>'+old_metafile:
                    result_value['metafile'] = '>'+new_metafile
                    if 'num_items' in result_value:
                        result_value['num_items'] = new_num_items
                else:
                    for additional_file, additional_output_file in zip(self.additional_dataset_files, self.additional_dataset_output_paths): #other .cs files in .csg
                        if result_value['metafile'] == '>'+additional_file:
                            new_additional_metafile = os.path.basename(additional_output_file)
                            result_value['metafile'] = '>' + new_additional_metafile
                            if 'num_items' in result_value:
                                result_value['num_items'] = new_num_items
        return group_data

    def write_csg(self, group_data, filepath):
        with open(filepath, 'w') as file:
            yaml.dump(group_data, file, default_flow_style=False)


    def open_filepicker(self, event, file_picker_queue):
        if self.relion_support and not self.cryosparc_support:
            title = "Select a relion star file"
            filetypes = [("Relion star file", "*.star"), ("All files", "*.*")]
        elif self.cryosparc_support and not self.relion_support:
            title = "Select a cryosparc csg/cs file"
            filetypes = [("Cryosparc group file", "*.csg"), ("Cryosparc data file", "*.cs"), ("All files", "*.*")]
        else:
            title = "Select a cryosparc csg/cs file or a relion star file"
            filetypes = [(("Cryosparc or Relion file", ("*.csg", "*.star"))), ("Cryosparc data file", "*.cs"), ("All files", "*.*")]
        file_path = filedialog.askopenfilename(title=title, filetypes=filetypes)
        if not (file_path == '' or file_path == () or file_path == [] or file_path == None):
            print_text("Selected file: %s" % file_path)
            self.dataset_path = file_path
            p.o_button.set_active(False)
            p.o_button_ax.set_visible(False)
            file_picker_queue.put(self.dataset_path)

    def get_poses(self):
        if self.is_cryosparc:
            poses = np.array(self.particles['alignments3D/pose'])
        elif self.is_relion:
            poses = np.array(self.particles[['rlnAngleRot', 'rlnAngleTilt', 'rlnAnglePsi']])
        return poses

    def mask_dataset(self, mask):
        mask_len = len(mask)
        if mask_len != self.num_particles:
            error_msg = 'Cannot apply mask of length %d to dataset of size %d' % (mask_len, self.num_particles)
            raise ValueError(error_msg)
        if self.is_cryosparc:
            return self.mask_cs(self.particles, mask)
        elif self.is_relion:
            return self.particles[mask]

    def mask_cs(self, particles, mask):
        return particles.mask(mask)

    def write_dataset(self, particles, file_path):
        if self.is_cryosparc:
            output_file_paths = self.write_cs_dataset(particles, file_path)
        elif self.is_relion:
            output_file_paths = self.write_rln_dataset(particles, file_path)
        else:
            output_file_paths = []
        return output_file_paths

    def write_cs_dataset(self, particles, file_path):
        output_csg = False if self.group_path == None else True
        output_file_paths = []
        num_particles = len(particles)
        if output_csg:
            group_data = self.update_csg(os.path.basename(file_path), num_particles)
            csg_filepath = file_path + 'g'
            self.write_csg(group_data, csg_filepath)
            output_file_paths.append(csg_filepath)
        self.write_cs(particles, file_path)
        output_file_paths.append(file_path)
        return output_file_paths

    def write_cs(self, particles, file_path):
        particles.save(file_path)

    def write_rln_dataset(self, particles, file_path):
        StarWriter.backup_if_file_exists = lambda *args: None #avoids a backup file error
        if self.optics is not None:
            starfile.write({'optics': self.optics, 'particles': particles}, file_path)
        else:
            starfile.write(particles, file_path)
        return [file_path]




class PlotObjects:
    def __init__(self, fig, root, canvas, plt, ax):
        self.fig = fig
        self.root = root
        self.canvas = canvas
        self.plt = plt
        self.ax = ax
        self.max_groups = 10

class CustomEllipse(patches.PathPatch):
    def __init__(self, center, width, height, angle=0.0, rotation_matrix_when_drawn=None, **kwargs):

        self.center = center
        self.width = width
        self.height = height
        self.angle = angle
        self.rotation_matrix_when_drawn = rotation_matrix_when_drawn
        self.y_upper_lim = np.pi / 2
        self.y_lower_lim = -np.pi / 2
        self.x_max_width = np.pi*2
        self.additional_lines = [None]

        ellipse_path = Path.unit_circle()
        transformation = Affine2D().scale(self.width / 2, self.height / 2).rotate_deg(self.angle).translate(
            self.center[0], self.center[1])
        ellipse_path = ellipse_path.transformed(transformation)

        #Cap path to not exceed y axis boundaries
        ellipse_path.vertices[:,1] = np.clip(ellipse_path.vertices[:,1],self.y_lower_lim, self.y_upper_lim)
        #prevent overlapping ellipses wider than 2 pi
        if self.x_max_width != None:
            ellipse_path.vertices[:,0] = np.clip(ellipse_path.vertices[:, 0], self.center[0]-(self.x_max_width/2.0), self.center[0]+(self.x_max_width/2.0))

        self.orig_ellipse_path = ellipse_path

        # Create a Path with the transformed ellipse path
        super().__init__(self.orig_ellipse_path, **kwargs)



class CustomTransform(Transform):
    def __init__(self, orig_rotation_matrix, rotation_matrix):
        Transform.__init__(self)
        self.orig_rotation_matrix = orig_rotation_matrix
        self.rotation_matrix = rotation_matrix

    def transform_non_affine(self, coords):
        #az = coords[:,0], el = coords[:,1]
        x = np.cos(coords[:,0]) * np.cos(coords[:,1])
        y = np.sin(coords[:,0]) * np.cos(coords[:,1])
        z = np.sin(coords[:,1])
        viewdirs = np.column_stack((x, y, z))
        reverse_matrix = np.matmul(np.linalg.inv(self.orig_rotation_matrix), self.rotation_matrix)
        viewdirs = np.matmul(viewdirs, reverse_matrix)
        az = np.arctan2(viewdirs[:, 1], viewdirs[:, 0])
        az = np.unwrap(az, discont=np.pi)
        el = np.arcsin(viewdirs[:, 2])
        #el = np.unwrap(el, discont=np.pi)
        coords = np.column_stack((az,el))
        return coords


class OffsetEllipseSelector(EllipseSelector):
    def __init__(self, offset, ax, onselect, callable_visual_props, **kwargs):
            self.offset = offset
            self.offset_ellipse1 = None
            self.offset_ellipse2 = None
            self.offset_ellipse1_extents = None
            self.offset_ellipse2_extents = None
            self.onselect_func = onselect
            self.callable_visual_props = callable_visual_props
            props = self.get_visual_props()
            super().__init__(ax, self.onselect, props=props, **kwargs)


    def onselect(self, eclick, erelease):
        self.onselect_func(eclick, erelease)
        if self.offset_ellipse1 in p.ax.patches:
            self.offset_ellipse1.remove()
        if self.offset_ellipse2 in p.ax.patches:
            self.offset_ellipse2.remove()
        #super().onselect(eclick, erelease)
        offset = self.offset

        x1, x2, y1, y2 = self.extents
        self.offset_ellipse1_extents = (x1, x2, y1, y2)
        self.offset_ellipse2_extents = (x1, x2, y1, y2)
        if supports_rotation:
            angle = self.rotation
        else:
            angle = 0.0
        if not ((x1 == x2) and (y1 == y2)):
            # Create and draw an offset ellipse
            self.offset_ellipse1_extents = (x1 - offset, x2 - offset, y1, y2)
            self.offset_ellipse2_extents = (x1 + offset, x2 + offset, y1, y2)
            props = self.get_visual_props()
            x1, x2, y1, y2 = self.offset_ellipse1_extents
            self.offset_ellipse1 = patches.Ellipse(((x1 + (x2 - x1) / 2.), y1 + (y2 - y1) / 2.), x2 - x1, y2 - y1, angle=angle, **props)
            x1, x2, y1, y2 = self.offset_ellipse2_extents
            self.offset_ellipse2 = patches.Ellipse(((x1 + (x2 - x1) / 2.), y1 + (y2 - y1) / 2.), x2 - x1, y2 - y1, angle=angle, **props)
            self.ax.add_patch(self.offset_ellipse1)
            self.ax.add_patch(self.offset_ellipse2)
            p.canvas.draw()

    def clear_all(self):
        self.extents = (0.0,0.0,0.0,0.0)
        self.clear()
        if self.offset_ellipse1 in p.ax.patches:
            self.offset_ellipse1.remove()
        if self.offset_ellipse2 in p.ax.patches:
            self.offset_ellipse2.remove()
        self.offset_ellipse1_extents = None
        self.offset_ellipse2_extents = None
        self.offset_ellipse1 = None
        self.offset_ellipse2 = None

    def get_visual_props(self):
        return {key: value() if callable(value) else value for key, value in self.callable_visual_props.items()}


def expmap(e, out=None):
    e = np.atleast_2d(e)
    theta = np.linalg.norm(e, axis=1)

    if out is None:
        out = np.eye(3, dtype=e.dtype).reshape(1, 3, 3).repeat(len(e), axis=0)

    zero_mask = theta < 1e-16
    nonzero_mask = ~zero_mask

    w = np.divide(e, theta[:, None], where=~zero_mask[:, None])

    s = np.sin(theta[nonzero_mask])
    c = 1 - np.cos(theta[nonzero_mask])

    w0 = s * w[nonzero_mask, 0]
    w1 = s * w[nonzero_mask, 1]
    w2 = s * w[nonzero_mask, 2]

    w00 = -w[nonzero_mask, 0] ** 2
    w11 = -w[nonzero_mask, 1] ** 2
    w22 = -w[nonzero_mask, 2] ** 2
    w01 = c * w[nonzero_mask, 0] * w[nonzero_mask, 1]
    w02 = c * w[nonzero_mask, 0] * w[nonzero_mask, 2]
    w12 = c * w[nonzero_mask, 1] * w[nonzero_mask, 2]

    out[nonzero_mask, 0, 0] = 1 + c * (w22 + w11)
    out[nonzero_mask, 0, 1] = w2 + w01
    out[nonzero_mask, 0, 2] = -w1 + w02
    out[nonzero_mask, 1, 0] = -w2 + w01
    out[nonzero_mask, 1, 1] = 1 + c * (w22 + w00)
    out[nonzero_mask, 1, 2] = w0 + w12
    out[nonzero_mask, 2, 0] = w1 + w02
    out[nonzero_mask, 2, 1] = -w0 + w12
    out[nonzero_mask, 2, 2] = 1 + c * (w11 + w00)

    # Handle cases where theta is close to zero
    out[zero_mask, 0, 0] = out[zero_mask, 1, 1] = out[zero_mask, 2, 2] = 1
    out[zero_mask, 0, 1:] = 0
    out[zero_mask, 1, 0] = out[zero_mask, 1, 2] = 0
    out[zero_mask, 2, :2] = 0

    return out


def matrix_from_euler(eulers, radians=False): # Create a rotation matrix from three Euler angles in ZYZ convention
    # eulers = [[rot, tilt, psi]]
    if eulers is not isinstance(eulers, np.ndarray):
        eulers = np.array(eulers)
    if eulers.ndim == 1:
        eulers = np.expand_dims(eulers, 0)
        one = True
    else:
        one = False
    num_items = eulers.shape[0]
    if radians == False:
        eulers = np.radians(eulers)
    m = np.zeros((num_items,3,3))
    m[:,0,0] = np.cos(eulers[:,2]) * np.cos(eulers[:,1]) * np.cos(eulers[:,0]) - np.sin(eulers[:,2]) * np.sin(eulers[:,0])
    m[:,0,1] = np.cos(eulers[:,2]) * np.cos(eulers[:,1]) * np.sin(eulers[:,0]) + np.sin(eulers[:,2]) * np.cos(eulers[:,0])
    m[:,0,2] = -np.cos(eulers[:,2]) * np.sin(eulers[:,1])
    m[:,1,0] = -np.sin(eulers[:,2]) * np.cos(eulers[:,1]) * np.cos(eulers[:,0]) - np.cos(eulers[:,2]) * np.sin(eulers[:,0])
    m[:,1,1] = -np.sin(eulers[:,2]) * np.cos(eulers[:,1]) * np.sin(eulers[:,0]) + np.cos(eulers[:,2]) * np.cos(eulers[:,0])
    m[:,1,2] = np.sin(eulers[:,2]) * np.sin(eulers[:,1])
    m[:,2,0] = np.sin(eulers[:,1]) * np.cos(eulers[:,0])
    m[:,2,1] = np.sin(eulers[:,1]) * np.sin(eulers[:,0])
    m[:,2,2] = np.cos(eulers[:,1])
    if one:
        m = np.squeeze(m, axis=0)
    return m

def euler_from_matrix(matrix, radians=False): #converts a matrix to Eulers
    if matrix is not isinstance(matrix, np.ndarray):
        matrix = np.array(matrix)
    if matrix.ndim == 2:
        matrix = np.expand_dims(matrix, 0)
        one = True
    else:
        one = False
    num_items = matrix.shape[0]
    eulers = np.zeros((num_items, 3))

    #tilt (0 to 180)
    eulers[:,1] = np.arccos(np.clip(matrix[:,2,2],-1.0,1.0))
    zero_mask = eulers[:,1] < 1e-16 #0.00001
    nonzero_mask = ~zero_mask

    #rot (-180 to 180)
    eulers[zero_mask, 0] = 0.0
    eulers[nonzero_mask, 0] = np.arctan2(matrix[nonzero_mask,2,1], matrix[nonzero_mask,2,0])

    #psi (-180 to 180)
    eulers[zero_mask, 2] = np.arctan2(-matrix[zero_mask,1,0], matrix[zero_mask,0,0])
    eulers[nonzero_mask, 2] = np.arctan2(matrix[nonzero_mask,1,2], -matrix[nonzero_mask,0,2])

    if radians == False:
        eulers = np.degrees(eulers)
    if one:
        eulers = np.squeeze(eulers, axis=0)
    return eulers

def identity_matrix():
    return np.eye(3)


def get_viewdir(df, rotation_matrix=None):
    Rs = df.rotation_matrices
    if rotation_matrix is not None:
        Rs = np.matmul(Rs, rotation_matrix)
    viewdirs = Rs[:,2]
    return viewdirs

def get_azimuth_elevation(viewdirs):
    azimuth = np.arctan2(viewdirs[:, 1], viewdirs[:, 0])
    elevation = np.arcsin(viewdirs[:, 2])
    return azimuth, elevation




def key_press_event(event):
    if event.key == 'enter':
        add_selection(event)
    if event.key == 'left':
        if supports_rotation:
            rotate_ellipse_left()
        else:
            print_text('EllipseSelector rotation not supported in this version of matplotlib.')
    if event.key == 'right':
        if supports_rotation:
            rotate_ellipse_right()
        else:
            print_text('EllipseSelector rotation not supported in this version of matplotlib.')

def select_callback(eclick, erelease):
    p.as_button.set_active(True)

def rotate_ellipse_left(step=10):
    p.selector.rotation = p.selector.rotation + step
    p.selector.onselect(None, None)

def rotate_ellipse_right(step=10):
    p.selector.rotation = p.selector.rotation - step
    p.selector.onselect(None, None)

def add_selection(event):
    current_group = radiobutton_get_state()
    x1,x2,y1,y2 = p.selector.extents
    if supports_rotation:
        angle = p.selector.rotation
    else:
        angle = 0.0
    if not (x1 == 0.0 and x2 == 0.0 and y1 == 0.0 and y2 == 0.0):
        if not ((x1 == x2) and (y1 == y2)):
            ellipse_triplet = [None, None, None]
            ellipse_triplet_extents = [p.selector.offset_ellipse1_extents, (x1, x2, y1, y2), p.selector.offset_ellipse2_extents]
            #ell = patches.Ellipse((x1+(x2-x1)/2.,y1+(y2-y1)/2.),x2-x1,y2-y1,angle=angle, fill=None, edgecolor=p.group_cmap(current_group), linewidth=6)
            for i,(ellipse, extents) in enumerate(zip(ellipse_triplet, ellipse_triplet_extents)):
                #offset = (i-1)*p.selector.offset
                if extents is not None:
                    x1, x2, y1, y2 = extents
                    ell = CustomEllipse(((x1+(x2-x1)/2.),y1+(y2-y1)/2.),x2-x1,y2-y1, angle=angle, fill=None, edgecolor=p.group_cmap(current_group), linewidth=6)
                    ell.rotation_matrix_when_drawn = p.rotation_matrix
                    p.ax.add_patch(ell)
                    ellipse_triplet[i] = ell
            p.selector.extents = (0,0,0,0)
            p.selector.clear_all()
            p.ellipse_groups[current_group].append(ellipse_triplet)
            p.ellipse_log.append(current_group)
            p.ag_button.set_active(True)
            p.f_button.set_active(True)
            p.u_button.set_active(True)
            p.fig.canvas.draw()
            print_text('Ellipse saved.')
            if len(p.ellipse_log) == 1:
                print_text('Draw another or click "New group" to make a selection that will be saved separately.')


def add_selection_group(event):
    if p.radio_buttons_visible == len(p.radio_buttons.labels):
        p.ag_button.set_active(False)
        p.fig.canvas.draw()
        return
    if p.radio_buttons_visible == 0:
        p.radio_buttons.circles[0].set_visible(True)
        p.radio_buttons.circles[1].set_visible(True)
        p.radio_buttons.labels[0].set_visible(True)
        p.radio_buttons.labels[1].set_visible(True)
        p.radio_buttons_visible = 2
        p.radio_buttons.set_active(1)
    elif p.radio_buttons_visible < len(p.radio_buttons.labels):
        next_button_index = p.radio_buttons_visible
        p.radio_buttons.circles[next_button_index].set_visible(True)
        p.radio_buttons.labels[next_button_index].set_visible(True)
        p.radio_buttons_visible = next_button_index+1
        p.radio_buttons.set_active(next_button_index)
    current_group = radiobutton_get_state()
    p.ellipse_groups.append([])
    p.selector.onselect(None,None)
    p.fig.canvas.draw()
    print_text('New selection group (%d) added.' % (current_group+1))
    if p.radio_buttons_visible == 2:
        print_text('Now draw another ellipse and hit enter (or click save selection).')
        print_text('(Any overlapping particles will default to lower numbered groups).')


def radiobuttonclick(event):
    selected = radiobutton_get_state()
    if selected > p.radio_buttons_visible-1:
        p.radio_buttons.set_active(p.radio_buttons.last_selected)
    else:
        p.radio_buttons.last_selected = selected
    current_group = radiobutton_get_state()
    p.selector.set_props(facecolor=p.group_cmap(current_group), edgecolor=p.group_cmap(current_group))
    p.fig.canvas.draw()

def radiobutton_get_state():
    return p.radio_labels.index(p.radio_buttons.value_selected)


def finish(event):
    if all([group == [] for group in p.ellipse_groups]):
        print_text('No selections made.')
        if save_image:
            print_text('Just outputting image of the plot.')
            dataset_basename, dataset_ext = os.path.splitext(df.dataset_path)
            image_path = dataset_basename + blank_image_string + '.png'
            p.plt.savefig(image_path)
            print_text(image_path)
    else:
        print_text('Outputting data...')
        p.as_button.set_active(False)
        p.ag_button.set_active(False)
        p.f_button.set_active(False)
        p.u_button.set_active(False)
        [button.set_active(False) for button in p.euler_buttons]
        p.finished = True
        output_thread.start()


def output_df():
    try:
        current_group = radiobutton_get_state()
        group_selections = []
        output_paths = []
        output_particle_numbers = []
        group_nums = []
        if df.is_cryosparc and df.additional_dataset_files is not []: #other .cs files in the .csg
            do_additional_files = True
            additional_dataset_masks = []
            additional_output_path_lists = []
        else:
            do_additional_files = False
        dataset_basename, dataset_ext = os.path.splitext(df.dataset_path)
        p.radio_buttons_visible = 1 if p.radio_buttons_visible == 0 else p.radio_buttons_visible
        for group in range(0,p.radio_buttons_visible):
            selected = np.zeros(df.num_particles)
            for ellipse_triplet in p.ellipse_groups[group]:
                for ellipse in ellipse_triplet:
                    if ellipse is not None:
                        path = ellipse.get_path()
                        if path is not None:
                            result = path.contains_points(np.column_stack((p.az, p.el)))
                            selected = np.logical_or(selected, result)
                # catch ellipses that have been transformed inside out due to both poles being surrounded
                if are_points_clockwise(ellipse_triplet[1].get_path().vertices[:-1]):
                    selected = np.invert(selected)
            if np.all(selected == False):
                continue
            for group_selection in group_selections: #overlapping particles default to the earliest group drawn
                selected = np.multiply(selected, np.invert(group_selection))
            group_selections.append(selected)
            subset = df.mask_dataset(selected)
            append_string = path_append_string + str(group+1) + dataset_ext
            subset_path = dataset_basename + append_string
            num_particles = len(subset)
            if do_additional_files:
                df.additional_dataset_output_paths = []
                for additional_file in df.additional_dataset_files:
                    additional_output_file = os.path.splitext(additional_file)[0]+append_string
                    df.additional_dataset_output_paths.append(additional_output_file)
                additional_output_path_lists.append(copy.copy(df.additional_dataset_output_paths))
                additional_dataset_masks.append(selected)
            paths = df.write_dataset(subset, subset_path)
            output_paths.extend(paths)
            if do_additional_files:
                output_paths.extend(df.additional_dataset_output_paths)
            group_nums.append(group+1)
            output_particle_numbers.append(num_particles)
        if len(group_selections) == 1:
            remainder = np.invert(group_selections[0])
        else:
            remainder = np.invert(np.logical_or.reduce(group_selections))
        subset = df.mask_dataset(remainder)
        append_string = remainder_string + dataset_ext
        remainder_path = dataset_basename + append_string
        num_particles = len(subset)
        if do_additional_files:
            df.additional_dataset_output_paths = []
            for additional_file in df.additional_dataset_files:
                additional_output_file = os.path.splitext(additional_file)[0] + append_string
                df.additional_dataset_output_paths.append(additional_output_file)
            additional_output_path_lists.append(copy.copy(df.additional_dataset_output_paths))
            additional_dataset_masks.append(remainder)
        paths = df.write_dataset(subset, remainder_path)
        output_paths.extend(paths)
        if do_additional_files:
            output_paths.extend(df.additional_dataset_output_paths)
        output_particle_numbers.append(num_particles)
        group_nums.append('remainder')
        p.selector.set_visible(False)
        p.selector.disconnect_events()
        if save_image:
            image_path = dataset_basename + image_string + '.png'
            p.plt.savefig(image_path)
            output_paths.append(image_path)
        if do_additional_files:
            df.particles = []
            for i,additional_file in enumerate(df.additional_dataset_files):
                file_root = os.path.dirname(df.dataset_path)
                additional_particles = df.load_cs(os.path.join(file_root,additional_file))
                for additional_mask, additional_output_files in zip(additional_dataset_masks, additional_output_path_lists):
                    additional_output_file = additional_output_files[i]
                    particles = df.mask_cs(additional_particles, additional_mask)
                    df.write_cs(particles, additional_output_file)
                additional_particles = []
        particle_total = 0
        print_text('Particle numbers:')
        for i,(output_particle_number,group_num) in enumerate(zip(output_particle_numbers, group_nums)):
            particle_total += output_particle_number
            if i != len(output_particle_numbers)-1:
                print_text("Group %d: %d" % (group_num, output_particle_number))
            else:
                print_text("Remainder: %d" % output_particle_number)
        print_text("Total number of particles: %d" % particle_total)
        print_text('Output file paths:')
        for output_path in output_paths:
            print_text(output_path)
        output_queue.put('Finished!')
    except Exception as ex:
        print_text(str(type(ex).__name__)+': '+str(ex), color='red')
        raise ex


def are_points_inside_ellipse(x_coords, y_coords, center, width, height, angle):
    # Convert to the local coordinate system of the ellipse
    x_local = x_coords - center[0]
    y_local = y_coords - center[1]

    # Apply the rotation matrix
    x_rotated = x_local * np.cos(angle) - y_local * np.sin(angle)
    y_rotated = x_local * np.sin(angle) + y_local * np.cos(angle)

    # Check the ellipse equation
    ellipse_equation = (x_rotated**2 / (width/2)**2) + (y_rotated**2 / (height/2)**2)
    return ellipse_equation <= 1

def are_points_clockwise(coords):
    n = coords.shape[0]
    if n < 3:
        return None #not enough points
    signed_area = np.sum(np.diff(coords[:,0]) * (coords[:-1,1] + coords[1:, 1]))
    if signed_area < 0: #anti-clockwise
        return False
    elif signed_area > 0: #clockwise
        return True
    else:
        return None #colinear

def get_y_intersect(coord1, coord2, y_axis_position):
    #y axis position must be > x of coord1 and < x of coord2
    x1, y1 = coord1
    x2, y2 = coord2
    if x1 == x2 or y_axis_position < min(x1, x2) or y_axis_position> max(x1,x2):
        return None
    intersect = ((y2-y1)/(x2-x1))*(y_axis_position-x1) + y1
    return np.array([y_axis_position, intersect])


def find_indices_crossing_value(array, value):
    cross_indices = np.where(np.diff(np.sign(array-value)))[0]
    indices_above = []
    indices_below = []

    for cross_index in cross_indices:
        index = cross_index #(np.arange(cross_index, cross_index - len(array), -1) % len(array))[0]
        indices_below.append(index)
        indices_above.append(index+1 % len(array))

    return indices_below, indices_above


def print_text(text, terminal_only=False, color=None):
    print(text)
    if not terminal_only:
        text_box.configure(state=tk.NORMAL)
        if color != None:
            text_box.tag_config('color', foreground=color)
            text_box.insert(tk.END, str(text) + "\n", 'color')
        else:
            text_box.insert(tk.END, str(text) + "\n")
        text_box.configure(state=tk.DISABLED)
        text_box.yview(tk.END)  # Auto-scroll to the end


def undo(event):
    last_group = p.ellipse_log[-1]
    if len(p.ellipse_groups[last_group]) != 0:
        ellipse_triplet = p.ellipse_groups[last_group][-1]
        for ell in ellipse_triplet:
            if ell is not None:
                for line in ell.additional_lines:
                    if line is not None:
                        if line in p.ax.lines:
                            line.remove()
                if ell in p.ax.patches:
                    ell.remove()
        del p.ellipse_groups[last_group][-1]
        del p.ellipse_log[-1]
    if len(p.ellipse_log) == 0:
        p.u_button.set_active(False)
    p.fig.canvas.draw()

def rotate_plot_button_event(event, euler_name, angle, queue, p):
    [button.set_active(False) for button in p.euler_buttons]
    as_button_status = p.as_button.get_active()
    ag_button_status = p.ag_button.get_active()
    u_button_status = p.u_button.get_active()
    f_button_status = p.f_button.get_active()
    p.as_button.set_active(False)
    p.ag_button.set_active(False)
    p.u_button.set_active(False)
    p.f_button.set_active(False)
    rot  = p.rot
    tilt = p.tilt
    psi = p.psi
    rot_min = -180.0
    rot_max = 180.0
    tilt_min = 0.0
    tilt_max = 180.0
    psi_min = -180.0
    psi_max = 180.0
    if euler_name == 'rot':
        rot += angle
        if rot <= rot_min:
            rot = rot_min
            print_text('Rot must be within %d and %d degrees.' % (rot_min, rot_max))
        elif rot >= rot_max:
            rot = rot_max
            print_text('Rot must be within %d and %d degrees.' % (rot_min, rot_max))
    elif euler_name == 'tilt':
        tilt += angle
        if tilt <= tilt_min:
            tilt = tilt_min
            print_text('Tilt must be within %d and %d degrees.' % (tilt_min, tilt_max))
        elif tilt >= tilt_max:
            tilt = tilt_max
            print_text('Tilt must be within %d and %d degrees.' % (tilt_min, tilt_max))
    elif euler_name == 'psi':
        psi += angle
        if psi <= psi_min:
            psi = psi_min
            print_text('Psi must be within %d and %d degrees.' % (psi_min, psi_max))
        elif psi >= psi_max:
            psi = psi_max
            print_text('Psi must be within %d and %d degrees.' % (psi_min, psi_max))
    else:
        raise ValueError('Unsupported euler_name %s' % str(euler_name))
    if not (rot == p.rot and tilt == p.tilt and psi == p.psi):
        p.rot = rot
        p.tilt = tilt
        p.psi = psi
        p.rot_text.set_text('Rot(z) %d°' % p.rot)
        p.tilt_text.set_text('Tilt(y) %d°' % p.tilt)
        p.psi_text.set_text('Psi(z) %d°' % p.psi)
        p.rotation_matrix = matrix_from_euler([rot, tilt, psi])
        print_text('Rotating plot...')
        queue.put('Rotate plot')
    else:
        p.rot = rot
        p.tilt = tilt
        p.psi = psi
        [button.set_active(True) for button in p.euler_buttons]
        p.as_button.set_active(as_button_status)
        p.ag_button.set_active(ag_button_status)
        p.u_button.set_active(u_button_status)
        p.f_button.set_active(f_button_status)

def rotate_plot(p, plot_data_queue):
    p.az, p.el = calculate_plot(df, p.rotation_matrix)
    plot_data_queue.put("Plot rotated")

def update_ellipse_plots(p): #update ellipse selections when the plot is rotated
    [patch.remove() for patch in p.ax.patches]
    [line.remove() for line in p.ax.lines]
    current_rotation_matrix = p.rotation_matrix
    for ellipse_group in p.ellipse_groups:
        for ellipse_triplet in ellipse_group:
            ellipse_paths = [None, None, None]
            ellipse = ellipse_triplet[1]
            ell_rotation_matrix = ellipse.rotation_matrix_when_drawn
            ellipse.set_path(ellipse.orig_ellipse_path)
            path = ellipse.get_path()
            custom_transform_obj = CustomTransform(ell_rotation_matrix, current_rotation_matrix)
            new_path = custom_transform_obj.transform_path(path)
            #catch the cases where an ellipse encloses a pole of the globe
            pole = 'none'  # north and south pole not in ellipse
            x_dif = new_path.vertices[-1, 0] - new_path.vertices[0, 0]
            if abs(x_dif) > np.pi:
                new_path = Path(np.flip(new_path.vertices[:-1], axis=0), closed=False)
                if x_dif >= 0:
                    pole = 'north'
                if x_dif < 0:
                    pole = 'south'
            if pole == 'none' and are_points_clockwise(new_path.vertices[:-1]):
                pole = 'both'
            ellipse_paths[1] = new_path
            for i, offset_ellipse in enumerate([ellipse_triplet[0], ellipse_triplet[2]]):
                if offset_ellipse is not None:
                    offset = -p.selector.offset if i ==0 else p.selector.offset
                    offset_path = new_path.transformed(Affine2D().translate(offset,0))
                    ellipse_paths[i*2] = offset_path
                else:
                    ellipse_paths[i * 2] = None
            ellipse_paths_without_none = [item for item in ellipse_paths if item is not None]
            if pole == 'north' or pole == 'south':
                if pole == 'north':
                    new_path = Path.make_compound_path(*ellipse_paths_without_none[::-1])
                    start_index = np.argmax(new_path.vertices[:, 0] < np.pi) - 1
                    end_index = np.argmax(new_path.vertices[:, 0] < -np.pi)
                    new_vertices = new_path.vertices[start_index:end_index + 1]
                    left_point = get_y_intersect(new_vertices[-2], new_vertices[-1], -np.pi)
                    right_point = get_y_intersect(new_vertices[0], new_vertices[1], np.pi)
                    new_vertices[-1] = left_point
                    new_vertices[0] = right_point
                    corner1 = [-np.pi, np.pi / 2.0]
                    corner2 = [np.pi, np.pi / 2.0]
                    new_vertices = np.vstack((corner2, new_vertices, corner1, corner2))
                if pole == 'south':
                    new_path = Path.make_compound_path(*ellipse_paths_without_none)
                    start_index = np.argmax(new_path.vertices[:, 0] > -np.pi) - 1
                    end_index = np.argmax(new_path.vertices[:, 0] > np.pi)
                    new_vertices = new_path.vertices[start_index:end_index + 1]
                    left_point = get_y_intersect(new_vertices[0], new_vertices[1], -np.pi)
                    right_point = get_y_intersect(new_vertices[-2], new_vertices[-1], np.pi)
                    new_vertices[0] = left_point
                    new_vertices[-1] = right_point
                    corner1 = [-np.pi, -np.pi / 2.0]
                    corner2 = [np.pi, -np.pi / 2.0]
                    new_vertices = np.vstack((corner1, new_vertices, corner2, corner1))
                new_path = Path(new_vertices,closed=True)
                ellipse_paths[1] = new_path
                ellipse_paths[0] = None
                ellipse_paths[2] = None
            elif pole == 'both': #add extra path around border
                left_upper_point = [-np.pi, 0]
                left_lower_point = [-np.pi, 0]
                right_upper_point = [np.pi, 0]
                right_lower_point = [np.pi, 0]
                for i,ellipse_path in enumerate(ellipse_paths):
                    if np.min(ellipse_path.vertices[:,0]) < -np.pi and np.max(ellipse_path.vertices[:,0]) > -np.pi: #is ellipse on left border
                        points = ellipse_path.vertices[:-2,:] #last two are redundant with first
                        indices_below, indices_above = find_indices_crossing_value(points[:,0], -np.pi)
                        if len(indices_below) == 1:
                            indices_above.append(indices_above[0])
                            indices_below.append(indices_below[0])
                        left_first_point = get_y_intersect(points[indices_below[0]], points[indices_above[0]],-np.pi)
                        left_second_point = get_y_intersect(points[indices_below[1]], points[indices_above[1]], -np.pi)
                        if left_first_point[1] < left_second_point[1]:
                            left_upper_point = left_second_point
                            left_lower_point = left_first_point
                        else:
                            left_upper_point = left_first_point
                            left_lower_point = left_second_point
                    if np.min(ellipse_path.vertices[:, 0]) < np.pi and np.max(ellipse_path.vertices[:, 0]) > np.pi:  # is ellipse on right border
                        points = ellipse_path.vertices[:-2, :]  # last two are redundant with first
                        indices_below, indices_above = find_indices_crossing_value(points[:,0], np.pi)
                        if len(indices_below) == 1:
                            indices_above.append(indices_above[0])
                            indices_below.append(indices_below[0])
                        right_first_point = get_y_intersect(points[indices_below[0]], points[indices_above[0]], np.pi)
                        right_second_point = get_y_intersect(points[indices_below[1]], points[indices_above[1]], np.pi)
                        if right_first_point[1] < right_second_point[1]:
                            right_upper_point = right_second_point
                            right_lower_point = right_first_point
                        else:
                            right_upper_point = right_first_point
                            right_lower_point = right_second_point
                upper_corner1 = [-np.pi, np.pi / 2.0]
                upper_corner2 = [np.pi, np.pi / 2.0]
                upper_vertices = np.vstack((left_upper_point, upper_corner1, upper_corner2, right_upper_point))
                lower_corner1 = [-np.pi, -np.pi / 2.0]
                lower_corner2 = [np.pi, -np.pi / 2.0]
                lower_vertices = np.vstack((left_lower_point, lower_corner1, lower_corner2, right_lower_point))
                ell = ellipse_triplet[1]
                upper_line = Line2D(*zip(*upper_vertices), color=ell.get_edgecolor(), linewidth=ell.get_linewidth())
                lower_line = Line2D(*zip(*lower_vertices), color=ell.get_edgecolor(), linewidth=ell.get_linewidth())
                p.ax.add_line(upper_line)
                p.ax.add_line(lower_line)
                ell.additional_lines = [upper_line, lower_line]
            for ellipse, ellipse_path in zip(ellipse_triplet, ellipse_paths):
                if ellipse is not None:
                    ellipse.set_path(ellipse_path)
                    if ellipse_path is not None:
                        p.ax.add_patch(ellipse)

def calculate_plot(df, rotation_matrix=None):
    if rotation_matrix is not None:
        p.rotation_matrix = rotation_matrix
    viewdirs = get_viewdir(df, rotation_matrix)
    az, el = get_azimuth_elevation(viewdirs)
    return az,el

def precalculate_plot(df):
    poses = df.get_poses()
    if df.is_cryosparc:
        df.rotation_matrices = expmap(poses)
    else:
        df.rotation_matrices = matrix_from_euler(poses, radians=False)


def generate_plot(df, p, plot_data_queue, file_picker_queue):
    try:
        if df.dataset_path == None:
            print_text('Open input file. (you can also pass a file directly on the command line).')
        while not p.stop_threads:
            try:
                queue_readout = file_picker_queue.get_nowait()
                if queue_readout == 'Rotate plot':
                    rotate_plot(p, plot_data_queue)
                else:
                    df.dataset_path = queue_readout
                    particles = df.load_dataset()
                    print_text('Calculating plot...')
                    precalculate_plot(df)
                    p.az,p.el = calculate_plot(df)
                    plot_data_queue.put("Plot ready")
            except queue.Empty:
                time.sleep(0.1)
    except Exception as ex:
        print_text(str(type(ex).__name__)+': '+str(ex), color='red')
        raise ex


def plot_data(p, plot_data_queue):
    try:
        queue_readout = plot_data_queue.get_nowait()
        if queue_readout == "Plot ready" or queue_readout == 'Plot rotated':
            print_text('Plotting...') if queue_readout == 'Plot ready' else None
            old_hb = p.hb if hasattr(p,'hb') else None
            p.hb = p.ax.hexbin(p.az, p.el, gridsize=50, bins='log', cmap=p.plot_cmap)
            if hasattr(p,'cb'):
                p.cb.remove()
            old_hb.remove() if old_hb != None else None
            p.cb = p.plt.colorbar(p.hb, ax=p.ax, label='# of images')
            update_ellipse_plots(p)
            if len(p.ellipse_log) > 0 and not p.finished:
                p.u_button.set_active(True)
                p.ag_button.set_active(True)
            if not p.finished:
                p.f_button.set_active(True)
            [button.set_active(True) for button in p.euler_buttons]
            print_text('Done.') if queue_readout == 'Plot ready' else None
            print_text('Draw an ellipse and hit enter (or click save selection).') if queue_readout == 'Plot ready' else None
            p.selector.clear_all()
            if supports_rotation:
                print_text('(You can rotate with left and right arrow keys).') if queue_readout == 'Plot ready' else None
            print_text('Plot rotated.') if queue_readout == 'Plot rotated' else None
            #simuluate_ellipses_test() if old_hb is None else None # uncomment for testing
            p.canvas.draw()
            root.after(100, lambda: plot_data(p, plot_data_queue))
    except queue.Empty:
        root.after(100, lambda: plot_data(p, plot_data_queue))


def output_finished(output_queue):
    try:
        queue_readout = output_queue.get_nowait()
        if queue_readout == 'Finished!':
            print_text('Finished!')
            output_thread.join(timeout=0.3)
            [button.set_active(True) for button in p.euler_buttons]
    except queue.Empty:
        root.after(100, lambda: output_finished(output_queue))


def on_close_window():
    print_text('Quitting...')
    p.stop_threads = True
    if thread.is_alive():
        thread.join(timeout=2)
    if output_thread.is_alive():
        output_thread.join(timeout=2)
    #root.destroy()
    root.quit()

def simuluate_ellipses_test(): #for testing only
    ellipses = [(-3.00984, -2.432200, 0.094789, 0.82602, 0.0),
    (-0.90194, -0.36483, 0.13541, 0.77185, 0.0),
    (1.155295, 1.78361, 0.12187, 0.8260, 0.0),
    (-2.01669, -1.327, -0.78539, -0.05416,0.0),
    (0.06080, 0.770196, -0.8260, 0.02708, 0.0),
    (2.2295, 2.8274, -0.7989, -0.05416, 0.0)]
    tilt_ellipses = [(0.4, -0.6, 0.4, -0.6, 20.0),
                     (2.8, 3.4, 0.4, -0.6, 0.0),
                    (-np.pi-0.2, 0.5, 0.3, -0.5, 0.0)]
    for ellipse in ellipses:
        p.selector.extents=ellipse[0:4]
        p.selector.rotation = ellipse[4]
        p.selector.onselect(None, None)
        time.sleep(0.2)
        add_selection(None)
        time.sleep(0.2)
        add_selection_group(None)
    rotate_plot_button_event(None, 'tilt', 90, file_picker_queue, p)
    time.sleep(0.5)
    for ellipse in tilt_ellipses:
        p.selector.extents=ellipse[0:4]
        p.selector.rotation = ellipse[4]
        p.selector.onselect(None, None)
        time.sleep(0.2)
        add_selection(None)
        time.sleep(0.2)
        add_selection_group(None)
    p.selector.clear_all()
    p.canvas.draw()


if __name__ == "__main__":
    argparser = ArgumentParser()
    args = argparser.parser.parse_args()
    argparser.validate(args)


    # Default values
    dataset_path = args.input_dataset_path
    path_append_string = '_group_'
    remainder_string = '_remainder'
    image_string = '_plot'
    blank_image_string = '_blank_plot'
    save_image = not args.dont_save_image
    group_cmap_name = args.group_cmap #'tab10'
    plot_cmap_name = args.plot_cmap #'jet'
    group_cmap = cm.get_cmap(group_cmap_name)

    print('Starting GUI...')
    #Start GUI
    root = tk.Tk()
    root.title("View Select - Select particles from a viewing direction distribution plot.")
    root.protocol("WM_DELETE_WINDOW", on_close_window)
    #Prepare GUI

    fig, ax = plt.subplots(figsize=(10, 4))

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    plt.subplots_adjust(right=0.9)

    # Create a scrolled text box
    text_box = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=110, height=10)
    text_box.pack(pady=10)

    # Set the text box to read-only initially
    text_box.configure(state=tk.DISABLED)

    # Set up the plot

    ax.set(xlim=(-math.pi, math.pi), ylim=(-math.pi / 2, math.pi / 2))
    ax.set_title("Viewing direction distribution")
    pi = np.pi
    ax.set_xticks(np.arange(-pi, pi + pi / 4, step=(pi / 4)),
                       ['-π', '-3π/4', '-π/2', '-π/4', '0', 'π/4', 'π/2', '3π/4', 'π'])
    ax.set_yticks(np.arange(-pi / 2, pi / 2 + pi / 4, step=(pi / 4)), ['-π/2', '-π/4', '0', 'π/4', 'π/2'])
    ax.set_xlabel('Azimuth')
    ax.set_ylabel('Elevation')

    p = PlotObjects(fig, root, canvas, plt, ax)

    p.group_cmap = group_cmap
    p.rotation_matrix = identity_matrix()
    p.rot = 0.0
    p.tilt = 0.0
    p.psi = 0.0
    p.angle_step = 90
    p.plot_cmap = plot_cmap_name

    p.ellipse_groups = []
    p.ellipse_groups.append([])
    p.ellipse_log = []
    p.finished = False
    p.stop_threads = False

    # Buttons
    plt.subplots_adjust(bottom=0.3)
    plt.subplots_adjust(left=0.18)
    as_button_ax = fig.add_axes([0.18, 0.08, 0.13, 0.07])  # [left, bottom, width, height]
    p.as_button = plt.Button(as_button_ax, 'Save selection')
    p.as_button.on_clicked(add_selection)
    p.as_button.set_active(False)

    ag_button_ax = fig.add_axes([0.33, 0.08, 0.13, 0.07])  # [left, bottom, width, height]
    p.ag_button = plt.Button(ag_button_ax, 'New group')
    p.ag_button.on_clicked(add_selection_group)
    p.ag_button.set_active(False)

    u_button_ax = fig.add_axes([0.48, 0.08, 0.13, 0.07])  # [left, bottom, width, height]
    p.u_button = plt.Button(u_button_ax, 'Undo')
    p.u_button.on_clicked(undo)
    p.u_button.set_active(False)

    f_button_ax = fig.add_axes([0.63, 0.08, 0.13, 0.07])  # [left, bottom, width, height]
    p.f_button = plt.Button(f_button_ax, 'Finish')
    p.f_button.on_clicked(finish)
    p.f_button.set_active(False)

    github_text = plt.text(0.995, 0.01, "github.com/robertstass/ViewSelect", color='grey', fontsize='medium', horizontalalignment='right', verticalalignment='bottom', transform=fig.transFigure)

    rot_neg_button_ax = fig.add_axes([0.02, 0.70, 0.04, 0.07])  # [left, bottom, width, height]
    p.rot_neg_button = plt.Button(rot_neg_button_ax, '-rot')
    p.rot_neg_button.on_clicked(lambda event: rotate_plot_button_event(event, 'rot', -p.angle_step, file_picker_queue, p))
    p.rot_neg_button.set_active(False)

    rot_pos_button_ax = fig.add_axes([0.07, 0.70, 0.04, 0.07])  # [left, bottom, width, height]
    p.rot_pos_button = plt.Button(rot_pos_button_ax, '+rot')
    p.rot_pos_button.on_clicked(lambda event: rotate_plot_button_event(event, 'rot', p.angle_step, file_picker_queue, p))
    p.rot_pos_button.set_active(False)
    p.rot_text = plt.text(0.025, 0.82, 'Rot(z) %d°' % p.rot, horizontalalignment='left', verticalalignment='center', transform=fig.transFigure)

    tilt_neg_button_ax = fig.add_axes([0.02, 0.5, 0.04, 0.07])  # [left, bottom, width, height]
    p.tilt_neg_button = plt.Button(tilt_neg_button_ax, '-tilt')
    p.tilt_neg_button.on_clicked(lambda event: rotate_plot_button_event(event, 'tilt', -p.angle_step, file_picker_queue, p))
    p.tilt_neg_button.set_active(False)

    tilt_pos_button_ax = fig.add_axes([0.07, 0.5, 0.04, 0.07])  # [left, bottom, width, height]
    p.tilt_pos_button = plt.Button(tilt_pos_button_ax, '+tilt')
    p.tilt_pos_button.on_clicked(lambda event: rotate_plot_button_event(event, 'tilt', p.angle_step, file_picker_queue, p))
    p.tilt_pos_button.set_active(False)
    p.tilt_text = plt.text(0.025, 0.62, 'Tilt(y) %d°' % p.tilt, horizontalalignment='left', verticalalignment='center', transform=fig.transFigure)

    psi_neg_button_ax = fig.add_axes([0.02, 0.3, 0.04, 0.07])  # [left, bottom, width, height]
    p.psi_neg_button = plt.Button(psi_neg_button_ax, '-psi')
    p.psi_neg_button.on_clicked(lambda event: rotate_plot_button_event(event, 'psi', -p.angle_step, file_picker_queue, p))
    p.psi_neg_button.set_active(False)

    psi_pos_button_ax = fig.add_axes([0.07, 0.3, 0.04, 0.07])  # [left, bottom, width, height]
    p.psi_pos_button = plt.Button(psi_pos_button_ax, '+psi')
    p.psi_pos_button.on_clicked(lambda event: rotate_plot_button_event(event, 'psi', p.angle_step, file_picker_queue, p))
    p.psi_pos_button.set_active(False)
    p.psi_text = plt.text(0.025, 0.42, 'Psi(z) %d°' % p.psi, horizontalalignment='left', verticalalignment='center', transform=fig.transFigure)

    p.euler_buttons = [p.rot_neg_button, p.rot_pos_button, p.tilt_neg_button, p.tilt_pos_button, p.psi_neg_button, p.psi_pos_button]
    


    # Radio button labels
    p.radio_labels = ['Group %d' % (i+1) for i in range(0,p.max_groups)]
    p.axRadio = plt.axes([0.82, 0.05, 0.4, 0.9])
    p.axRadio.axis('off')
    p.radio_buttons = RadioButtons(p.axRadio, p.radio_labels, activecolor='k')
    p.radio_buttons.on_clicked(radiobuttonclick)
    p.radio_buttons_visible = 0
    p.radio_buttons.last_selected = 0
    for i, (circle, label) in enumerate(zip(p.radio_buttons.circles, p.radio_buttons.labels)):
        circle.set_radius(0.035)
        circle.set_edgecolor(group_cmap(i))
        circle.set_linewidth(2)
        circle.set_visible(False)
        label.set_visible(False)


    p.selector_props = {'facecolor': lambda : group_cmap(radiobutton_get_state()), 'edgecolor': lambda : group_cmap(radiobutton_get_state()), 'alpha': 0.3, 'fill': True}
    #props = {key: value() if callable(value) else value for key, value in p.selector_props.items()}
    p.selector = OffsetEllipseSelector(np.pi * 2, ax, select_callback, callable_visual_props=p.selector_props, useblit=True, button=[1], minspanx=5, minspany=5,
                spancoords='pixels', ignore_event_outside=False, interactive=True)
    supports_rotation = True if hasattr(p.selector, 'rotation') else False


    # key bindings
    fig.canvas.mpl_connect('key_press_event', key_press_event)


    #Dataset
    df = CryoDatafile(dataset_path, cryosparc_support=cryosparc_support, relion_support=relion_support)


    #Queues

    plot_data_queue = queue.Queue()
    file_picker_queue = queue.Queue()
    output_queue = queue.Queue()

    #Open file picker

    if dataset_path == None:
        p.o_button_ax = fig.add_axes([0.35, 0.5, 0.3, 0.2])  # [left, bottom, width, height]
        p.o_button = plt.Button(p.o_button_ax, 'Open')
        p.o_button.on_clicked(lambda event: df.open_filepicker(event, file_picker_queue))
    else:
        file_picker_queue.put(dataset_path)


    #Threading
    try:
        stop_event = Event()
        thread = Thread(target=generate_plot, args=(df, p, plot_data_queue, file_picker_queue))
        output_thread = Thread(target=output_df)
        thread.start()

        root.after(500, lambda: plot_data(p, plot_data_queue))
        root.after(500, lambda: output_finished(output_queue))


        root.mainloop()

    except Exception as ex:
        print_text(str(type(ex).__name__)+': '+str(ex), color='red')
        raise ex

    stop_event.set()
    p.stop_threads = True
    if thread.is_alive():
        thread.join(timeout=0.1)
    if output_thread.is_alive():
        output_thread.join(timeout=0.1)




