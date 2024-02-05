#!/usr/bin/env python
import os
import sys
import math
import numpy as np
from threading import Thread, Event
import argparse
import queue
import copy

import matplotlib as mpl
from cycler import cycler
from matplotlib import pyplot as plt
from matplotlib.widgets import  EllipseSelector, Slider, Button, RadioButtons, TextBox
from matplotlib.transforms import Transform
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.transforms import Affine2D

import tkinter as tk
from tkinter import scrolledtext
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog

from cryosparc.tools import CryoSPARC
from cryosparc.dataset import Dataset
import yaml

import warnings
warnings.filterwarnings("ignore", category=mpl.MatplotlibDeprecationWarning)
#mpl.use("Qt5Agg")




class ArgumentParser():

    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''
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
            "conda activate view_select" followed by "python view_select.py"
            
            Why do this?
            For small and/or featureless samples, it is common to get one or two views that dominate all others resulting
            in a preferred orientation problem. It is then optimal to try to enrich other views but these can sometimes be
            hard to find at the 2D classification stage. The 3D refinement stage can sometimes help but it can be difficult
            to relate hotspots in a viewing direction distribution plot back to the 2D classes. There is a partial solution
            to this which is to generate a .bild file for viewing in chimera using pyem (https://github.com/asarnow/pyem).
            However, "View Select" allows the user to select particles from the plot that can then be used for further 2D classification
            and/or particle picking with tools like Topaz.
            
            Written by Robert Stass, Bowden group, STRUBI/OPIC (2024)
            
            ''')
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        add('--dont_save_image', action='store_true', help='Stop the script outputting an image at the end.')
        add('input_dataset_path', nargs='?', default=None, help='File path of an exported cryosparc .csg or .cs file. Can leave blank to open a file in the GUI.')

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












######### Main script ##########


class CryosparcDataset:
    def __init__(self, dataset_path=None, group_path=None):
        self.group_path = group_path #.csg
        self.dataset_path = dataset_path #.cs
        self.num_particles = None

    def load_dataset(self):
        splitpath = os.path.splitext(self.dataset_path)
        if splitpath[1] == '.csg':
            print_text('.csg file supplied.')
            print_text('Reading .csg file...')
            self.group_path = self.dataset_path
            self.load_csg()
        else:
            print_text('.cs file supplied. (no .csg file will be produced).')
        print_text('Reading .cs file...')
        self.particles_dset = Dataset.load(self.dataset_path)
        self.num_particles = len(self.particles_dset)
        return self.particles_dset

    def load_csg(self):
        with open(self.group_path, 'r') as file:
            self.group_data = yaml.safe_load(file)
        cs = self.group_data.get('results', {}).get('alignments3D', {}).get('metafile')
        if cs == None:
            error_msg = 'Error: .csg file must contain an "alignments3D" field.'
            print_text(error_msg)
            raise OSError(error_msg)
        if cs is not None and cs.startswith('>'):
            cs = cs[1:]
        dirname = os.path.dirname(self.group_path)
        if dirname != '':
            cs = os.path.join(dirname, cs)
        if os.path.exists(cs):
            self.dataset_path = cs
        else:
            error_msg = 'Error: When supplying a .csg file, the corresponding .cs file (%s) must be present in the same directory.' % (cs)
            print_text(error_msg)
            raise OSError(error_msg)

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
        return group_data

    def write_csg(self, group_data, filepath):
        with open(filepath, 'w') as file:
            yaml.dump(group_data, file, default_flow_style=False)


    def open_filepicker(self, event, file_picker_queue):
        file_path = filedialog.askopenfilename(title="Select a cryosparc file",
                                               filetypes=[("Cryosparc group file", "*.csg"), ("Cryosparc data file", "*.cs"), ("All files", "*.*")])
        if not (file_path == '' or file_path == () or file_path == [] or file_path == None):
            print("Selected file:", file_path)
            self.dataset_path = file_path
            p.o_button.set_active(False)
            p.o_button_ax.set_visible(False)
            file_picker_queue.put(self.dataset_path)


class PlotObjects:
    def __init__(self, fig, root, canvas, plt, ax):
        self.fig = fig
        self.root = root
        self.canvas = canvas
        self.plt = plt
        self.ax = ax
        #self.current_group = 0
        self.max_groups = 10

class CustomEllipse(mpl.patches.PathPatch):
    def __init__(self, center, width, height, angle=0.0, rotation_matrix_when_drawn=None, **kwargs):

        self.center = center
        self.width = width
        self.height = height
        self.angle = angle
        self.rotation_matrix_when_drawn = rotation_matrix_when_drawn

        ellipse_path = Path.unit_circle()
        transformation = Affine2D().scale(self.width / 2, self.height / 2).rotate_deg(self.angle).translate(
            self.center[0], self.center[1])
        self.orig_ellipse_path = ellipse_path.transformed(transformation)

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
        #reverse_matrix = np.linalg.inv(self.rotation_matrix)
        reverse_matrix = np.matmul(self.orig_rotation_matrix, self.rotation_matrix)
        reverse_matrix = np.matmul(reverse_matrix, self.orig_rotation_matrix)
        viewdirs = np.matmul(viewdirs, reverse_matrix)
        az = np.arctan2(viewdirs[:, 1], viewdirs[:, 0])
        az = ((az+2*np.pi) % (2*np.pi))-(2*np.pi)
        el = np.arcsin(viewdirs[:, 2])
        #az, el = get_azimuth_elevation(viewdirs)
        coords = np.column_stack((az,el))
        return coords


class OffsetEllipseSelector(EllipseSelector):
    def __init__(self, offset, ax, onselect, **kwargs):
            self.offset = offset
            self.offset_ellipse1 = None
            self.offset_ellipse2 = None
            self.onselect_func = onselect
            super().__init__(ax, self.onselect, **kwargs)


    def onselect(self, eclick, erelease):
        self.onselect_func(eclick, erelease)
        if self.offset_ellipse1:
            self.offset_ellipse1.remove()
        if self.offset_ellipse2:
            self.offset_ellipse2.remove()
        #super().onselect(eclick, erelease)
        offset = self.offset

        x1, x2, y1, y2 = self.extents
        if supports_rotation:
            angle = self.rotation
        else:
            angle = 0.0
        if not (x1 == 0.0 and x2 == 0.0 and y1 == 0.0 and y2 == 0.0):
            if not ((x1 == x2) and (y1 == y2)):
                # Create and draw an offset ellipse
                self.offset_ellipse1 = mpl.patches.Ellipse(((x1 + (x2 - x1) / 2.) - offset, y1 + (y2 - y1) / 2.), x2 - x1, y2 - y1, angle=angle, **self._props)
                self.offset_ellipse2 = mpl.patches.Ellipse(((x1 + (x2 - x1) / 2.) + offset, y1 + (y2 - y1) / 2.), x2 - x1, y2 - y1, angle=angle, **self._props)
                self.ax.add_patch(self.offset_ellipse1)
                self.ax.add_patch(self.offset_ellipse2)
                #p.canvas.draw()
                #self.offset_ellipse1 = None  # Reset offset ellipse after selection
                #self.offset_ellipse1 = None
    def clear_all(self):
        self.clear()
        if self.offset_ellipse1:
            self.offset_ellipse1.remove()
        if self.offset_ellipse2:
            self.offset_ellipse2.remove()
        self.offset_ellipse1 = None
        self.offset_ellipse2 = None



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

def identity_matrix():
    return np.eye(3)

def Rx(theta):
    return np.array([[1, 0, 0],
                      [0, math.cos(theta), -math.sin(theta)],
                      [0, math.sin(theta), math.cos(theta)]])

def Ry(theta):
    return np.array([[math.cos(theta), 0, math.sin(theta)],
                      [0, 1, 0],
                      [-math.sin(theta), 0, math.cos(theta)]])

def Rz(theta):
    return np.array([[math.cos(theta), -math.sin(theta), 0],
                      [math.sin(theta), math.cos(theta), 0],
                      [0, 0, 1]])

def get_viewdir(poses, rotation_matrix=None):
    Rs = expmap(poses)
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
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    p.as_button.set_active(True)

def rotate_ellipse_left(step=10):
    p.selector.rotation = p.selector.rotation + step

def rotate_ellipse_right(step=10):
    p.selector.rotation = p.selector.rotation - step

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
            #ell = mpl.patches.Ellipse((x1+(x2-x1)/2.,y1+(y2-y1)/2.),x2-x1,y2-y1,angle=angle, fill=None, edgecolor=p.group_cmap(current_group), linewidth=6)
            for i,ellipse in enumerate(ellipse_triplet):
                offset = (i-1)*p.selector.offset
                ell = CustomEllipse(((x1+(x2-x1)/2.)+offset,y1+(y2-y1)/2.),x2-x1,y2-y1, angle=angle, fill=None, edgecolor=p.group_cmap(current_group), linewidth=6)
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
            dataset_basename, dataset_ext = os.path.splitext(cs.dataset_path)
            image_path = dataset_basename + blank_image_string + '.png'
            p.plt.savefig(image_path)
            print_text(image_path)
    else:
        print_text('Outputting data...')
        p.as_button.set_active(False)
        p.ag_button.set_active(False)
        p.f_button.set_active(False)
        p.u_button.set_active(False)
        output_thread.start()


def output_cs():
    try:
        current_group = radiobutton_get_state()
        group_selections = []
        output_paths = []
        output_particle_numbers = []
        group_nums = []
        output_csg = False if cs.group_path == None else True
        dataset_basename, dataset_ext = os.path.splitext(cs.dataset_path)
        p.radio_buttons_visible = 1 if p.radio_buttons_visible == 0 else p.radio_buttons_visible
        for group in range(0,p.radio_buttons_visible):
            selected = np.zeros(cs.num_particles)
            for ellipse in p.ellipse_groups[group]:
                xy = p.hb.get_offsets()
                center = ellipse.center
                result = are_points_inside_ellipse(p.az, p.el, ellipse.center, ellipse.width, ellipse.height, ellipse.angle)
                selected = np.logical_or(selected, result)
            if np.all(selected == False):
                continue
            for group_selection in group_selections: #overlapping particles default to the earliest group drawn
                selected = np.multiply(selected, np.invert(group_selection))
            group_selections.append(selected)
            subset = cs.particles_dset.mask(selected)
            subset_path = dataset_basename + path_append_string + str(group+1) + dataset_ext
            num_particles = len(subset)
            if output_csg:
                group_data = cs.update_csg(os.path.basename(subset_path), num_particles)
                csg_filepath = subset_path + 'g'
                cs.write_csg(group_data, csg_filepath)
                output_paths.append(csg_filepath)
            group_nums.append(group+1)
            output_particle_numbers.append(num_particles)
            output_paths.append(subset_path)
            subset.save(subset_path)
        if len(group_selections) == 1:
            remainder = np.invert(group_selections[0])
        else:
            remainder = np.invert(np.logical_or.reduce(group_selections))
        subset = cs.particles_dset.mask(remainder)
        remainder_path = dataset_basename + remainder_string + dataset_ext
        num_particles = len(subset)
        if output_csg:
            group_data = cs.update_csg(os.path.basename(remainder_path), num_particles)
            csg_filepath = remainder_path + 'g'
            cs.write_csg(group_data, csg_filepath)
            output_paths.append(csg_filepath)
        output_particle_numbers.append(num_particles)
        output_paths.append(remainder_path)
        group_nums.append('remainder')
        subset.save(remainder_path)
        p.selector.set_visible(False)
        p.selector.disconnect_events()
        if save_image:
            image_path = dataset_basename + image_string + '.png'
            p.plt.savefig(image_path)
            output_paths.append(image_path)
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
        #print_text('Finished!')
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

# Create a function to append text to the text box
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
        ell = p.ellipse_groups[last_group][-1]
        ell.remove()
        del p.ellipse_groups[last_group][-1]
        del p.ellipse_log[-1]
    if len(p.ellipse_log) == 0:
        p.u_button.set_active(False)
    p.fig.canvas.draw()

def rotate_plot(event, axis, angle, plot_data_queue, p):
    if axis == 'x':
        rotation_matrix = Rx(int(angle))
    elif axis == 'y':
        rotation_matrix = Ry(int(angle))
    elif axis == 'z':
        rotation_matrix = Rz(int(angle))
    else:
        raise ValueError('Unsupported axis %s' % str(axis))
    p.az, p.el = calculate_plot(cs.particles_dset, rotation_matrix)
    plot_data_queue.put("Plot ready")

def update_ellipse_plots(p):
    current_rotation_matrix = p.rotation_matrix
    for ellipse_group in p.ellipse_groups:
        for ellipse_triplet in ellipse_group:
            ellipse = ellipse_triplet[1]
            ell_rotation_matrix = ellipse.rotation_matrix_when_drawn
            path = ellipse.get_path()
            ellipse.set_path(ellipse.orig_ellipse_path)
            path = ellipse.get_path()
            print(ell_rotation_matrix)
            print(current_rotation_matrix)
            #transform_rotation_matrix = np.matmul(ell_rotation_matrix, current_rotation_matrix)
            #current_rotation_matrix = np.linalg.inv(current_rotation_matrix)
            custom_transform_obj = CustomTransform(ell_rotation_matrix, current_rotation_matrix)
            new_path = custom_transform_obj.transform_path(path)
            ellipse.set_path(new_path)
            ellipse.remove()
            p.ax.add_patch(ellipse)
            for i, offset_ellipse in enumerate([ellipse_triplet[0], ellipse_triplet[2]]):
                offset = -p.selector.offset if i ==0 else p.selector.offset
                offset_path = new_path.transformed(Affine2D().translate(offset,0))
                offset_ellipse.set_path(offset_path)
                offset_ellipse.remove()
                p.ax.add_patch(offset_ellipse)

def calculate_plot(particles_dset, rotation_matrix=None):
    print_text('Calculating plot...')
    poses = np.array(particles_dset['alignments3D/pose'])
    if rotation_matrix is not None:
        rotation_matrix = np.matmul(p.rotation_matrix, rotation_matrix)
        p.rotation_matrix = rotation_matrix
    viewdirs = get_viewdir(poses, rotation_matrix)
    az, el = get_azimuth_elevation(viewdirs)
    return az,el


def generate_plot(cs, p, plot_data_queue, file_picker_queue):
    try:
        if cs.dataset_path == None:
            print_text('Open input file. (you can also pass a file directly on the command line).')
        cs.dataset_path = file_picker_queue.get()
        particles_dset = cs.load_dataset()
        p.az,p.el = calculate_plot(particles_dset)
        plot_data_queue.put("Plot ready")
        return
    except Exception as ex:
        print_text(str(type(ex).__name__)+': '+str(ex), color='red')
        raise ex


def plot_data(p, plot_data_queue):
    try:
        queue_readout = plot_data_queue.get_nowait()
        if queue_readout == "Plot ready":
            print_text('Plotting...')
            old_hb = p.hb if hasattr(p,'hb') else None
            p.hb = p.ax.hexbin(p.az, p.el, gridsize=50, bins='log', cmap='jet')
            if hasattr(p,'cb'):
                p.cb.remove()
            old_hb.remove() if old_hb != None else None
            p.cb = p.plt.colorbar(p.hb, ax=p.ax, label='# of images')
            update_ellipse_plots(p)
            p.canvas.draw()
            p.f_button.set_active(True)
            p.x_neg_90_button.set_active(True)
            print_text('Done.')
            print_text('Draw an ellipse and hit enter (or click save selection).')
            if supports_rotation:
                print_text('(You can rotate with left and right arrow keys).')
            root.after(100, lambda: plot_data(p, plot_data_queue))
    except queue.Empty:
        root.after(100, lambda: plot_data(p, plot_data_queue))


def output_finished(output_queue):
    try:
        queue_readout = output_queue.get_nowait()
        if queue_readout == 'Finished!':
            print_text('Finished!')
    except queue.Empty:
        root.after(100, lambda: output_finished(output_queue))


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
    colormap_name = 'tab10'
    try:
        group_cmap = mpl.colormaps.get_cmap(colormap_name)
    except AttributeError:
        group_cmap = mpl.cm.get_cmap(colormap_name)
    
    #Start GUI
    root = tk.Tk()
    root.title("View Select - Select particles from a viewing direction distribution plot.")
    
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

    p.selector = OffsetEllipseSelector(np.pi*2, ax, select_callback, useblit=True, button=[1], minspanx=5, minspany=5, spancoords='pixels',
        props={'facecolor': group_cmap(0), 'edgecolor': group_cmap(0), 'alpha': 0.3, 'fill': True},
        ignore_event_outside=False, interactive=True)
    #p.selector.state_modifier_keys["rotate"] = 'alt'
    supports_rotation = True if hasattr(p.selector, 'rotation') else False
    #p.selector.

    p.ellipse_groups = []
    p.ellipse_groups.append([])
    p.ellipse_log = []

    # Buttons
    plt.subplots_adjust(bottom=0.3)
    as_button_ax = fig.add_axes([0.1, 0.05, 0.13, 0.05])  # [left, bottom, width, height]
    p.as_button = plt.Button(as_button_ax, 'Save selection')
    p.as_button.on_clicked(add_selection)
    p.as_button.set_active(False)

    ag_button_ax = fig.add_axes([0.25, 0.05, 0.13, 0.05])  # [left, bottom, width, height]
    p.ag_button = plt.Button(ag_button_ax, 'New group')
    p.ag_button.on_clicked(add_selection_group)
    p.ag_button.set_active(False)

    u_button_ax = fig.add_axes([0.40, 0.05, 0.13, 0.05])  # [left, bottom, width, height]
    p.u_button = plt.Button(u_button_ax, 'Undo')
    p.u_button.on_clicked(undo)
    p.u_button.set_active(False)

    f_button_ax = fig.add_axes([0.55, 0.05, 0.13, 0.05])  # [left, bottom, width, height]
    p.f_button = plt.Button(f_button_ax, 'Finish')
    p.f_button.on_clicked(finish)
    p.f_button.set_active(False)

    x_neg_90_button_ax = fig.add_axes([0.02, 0.8, 0.08, 0.05])  # [left, bottom, width, height]
    p.x_neg_90_button = plt.Button(x_neg_90_button_ax, 'Rot x -90')
    p.x_neg_90_button.on_clicked(lambda event: rotate_plot(event, 'x', -90, plot_data_queue, p))
    p.x_neg_90_button.set_active(False)

    


    # Radio button labels
    p.radio_labels = ['Group %d' % (i+1) for i in range(0,p.max_groups)]
    p.axRadio = plt.axes([0.815, 0.05, 0.4, 0.9])
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


    # key bindings
    fig.canvas.mpl_connect('key_press_event', key_press_event)


    #Dataset
    cs = CryosparcDataset(dataset_path)



    #Queues

    plot_data_queue = queue.Queue()
    file_picker_queue = queue.Queue()
    output_queue = queue.Queue()

    #Open file picker

    if dataset_path == None:
        p.o_button_ax = fig.add_axes([0.35, 0.5, 0.3, 0.2])  # [left, bottom, width, height]
        p.o_button = plt.Button(p.o_button_ax, 'Open')
        p.o_button.on_clicked(lambda event: cs.open_filepicker(event, file_picker_queue))
    else:
        file_picker_queue.put(dataset_path)


    #Threading
    try:
        stop_event = Event()
        thread = Thread(target=generate_plot, args=(cs, p, plot_data_queue, file_picker_queue))
        output_thread = Thread(target=output_cs)
        thread.start()

        root.after(500, lambda: plot_data(p, plot_data_queue))
        root.after(500, lambda: output_finished(output_queue))


        root.mainloop()

    except Exception as ex:
        print_text(str(type(ex).__name__)+': '+str(ex), color='red')
        raise ex

    stop_event.set()
    thread.join()
    output_thread.join()





