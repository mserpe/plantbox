import sys; sys.path.append("../src/python_modules/"); sys.path.append("../")

import os
import pandas as pd
import numpy as np
import argparse

import viewer_conductivities
from viewer_data import ViewerDataModel
import xylem_flux

"""
creates a csv file per rsml file containing krs values, and suf values per 1 mm layers

expected three arguments: file path, scenario index (1-4), (optionally, --shoot, only use for measured monocots)

e.g. 
python3 label_rsml.py ~/workspace/DUMUX/CPlantBox/gui/maize 1
python3 label_rsml.py ~/workspace/DUMUX/CPlantBox/gui/dicot 1

"""


def label_file(file, artificial_shoot, scenario_index):
    """ opens rsml and computes krs and suf """
    data = ViewerDataModel()
    data.open_rsml(file)  # same as in gui (to check)
    if artificial_shoot:
        data.add_artificial_shoot()

    r = data.xylem_flux
    j = scenario_index - 1
    if j == 0:
        viewer_conductivities.init_constant_scenario1(r)
    elif j == 1:
        viewer_conductivities.init_constant_scenario2(r)
    elif j == 2:
        viewer_conductivities.init_dynamic_scenario1(r)
    elif j == 3:
        viewer_conductivities.init_dynamic_scenario2(r)

    krs, _ = data.xylem_flux.get_krs(data.max_ct, data.base_segs)
    nop = len(data.base_nodes)  # number of plants (we might want to multiply suf by it?)
    nop = 1
    suf = r.get_suf(data.max_ct) * nop
    data.analyser.addData("SUF", suf)
    n = int(np.ceil(-data.analyser.getMinBounds().z))
    suf_ = data.analyser.distribution("SUF", 0., float(-n), int(n) * 10, False)  # mm layers
    depth = r.get_mean_suf_depth(data.max_ct)
    print("depth: ", depth, "cm")
    return suf_, krs, depth


def write_csv(file, suf_, krs, si, depth):
    """ writes an xls sheet containing suf """
    file_csv = file.rsplit('.', 1)[0] + '_suf' + str(si) + '.csv'
    print("krs {:g}, suf from {:g} to {:g}, sums up to {:g}".format(krs, np.min(suf_), np.max(suf_), np.sum(suf_)))
    suf_ = np.insert(suf_, 0, depth)
    suf_ = np.insert(suf_, 0, krs)
    df = pd.DataFrame(suf_)
    df.to_csv(file_csv, index = False, header = False)
    print("done. \n\n")


if __name__ == '__main__':

    """ parse arguments """
    parser = argparse.ArgumentParser()
    parser.add_argument('file_path', type = str, help = 'file path')
    parser.add_argument('scenario_index', type = int, help = 'scenario index (1-4)')
    parser.add_argument('--shoot', action = 'store_true', help = 'adds an artificial shoot')
    args = parser.parse_args()

    walk_dir = args.file_path
    print('walk_dir = ' + walk_dir)

    # If your current working directory may change during script execution, it's recommended to
    # immediately convert program arguments to an absolute path. Then the variable root below will
    # be an absolute path as well. Example:
    # walk_dir = os.path.abspath(walk_dir)
    print('walk_dir (absolute) = ' + os.path.abspath(walk_dir))
    artificial_shoot = args.shoot

    scenario_index = args.scenario_index
    print("Scenario index {:g}, see file viewer_conductivities.py".format(scenario_index))

    if scenario_index < 1 or scenario_index > 4:
        raise

    print()

    """ walk the files """
    for root, subdirs, files in os.walk(walk_dir):
        # print('--\nroot = ' + root)
        list_file_path = os.path.join(root, 'my-directory-list.txt')
        # print('list_file_path = ' + list_file_path)

        with open(list_file_path, 'wb') as list_file:

            for subdir in subdirs:
                print('\t- subdirectory ' + subdir)

            for filename in files:

                if filename.endswith('.rsml'):

                    file_path = os.path.join(root, filename)
                    print('file %s (full path: %s)\n' % (filename, file_path))

                    suf_, krs, depth = label_file(file_path, artificial_shoot, scenario_index)
                    write_csv(file_path, suf_, krs, scenario_index, depth)

