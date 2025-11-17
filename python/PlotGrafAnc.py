#!/usr/bin/env python3

import argparse

import os
import string
import sys

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.patches as patches

import matplotlib.colors as mcolors
from pathlib import Path
from matplotlib import colormaps
from matplotlib import ticker

from PopColors import PopColor
from GraphParameters import GraphParameter
from SubjectPops import SubjectPops
from matplotlib.patches import Ellipse

def ShowLegends(ax, graph_params, sbj_pop_obj):
    """
    Display legends on the user-specified position
    """
    x_range = graph_params.xmax - graph_params.xmin
    y_range = graph_params.ymax - graph_params.ymin

    lgd_width  = graph_params.lgd_size * x_range * 0.1 / graph_params.width
    lgd_height = graph_params.lgd_size * y_range * 0.1 / graph_params.height

    lgd_x_val = graph_params.xmin + graph_params.lgd_x_pos * x_range
    lgd_y_val = graph_params.ymax - graph_params.lgd_y_pos * y_range
    
    lgd_x_gap = lgd_width  * graph_params.lgd_x_gap
    lgd_y_gap = lgd_height * graph_params.lgd_y_gap

    num_lgds = len(sbj_pop_obj.lgd_pops)

    x = lgd_x_val
    row = 0

    for i in range(num_lgds):
        lbl = sbj_pop_obj.lgd_pops[i]
        if lbl == "":
            lbl = "NO VALUE"
            
        y = lgd_y_val - lgd_y_gap * row
        cx, cy = x + lgd_width/2, y + lgd_height/2
        radius = lgd_width/2
            
        color = sbj_pop_obj.lgd_colors[i]
        try:
            ellipse = Ellipse(xy=(cx, cy), width=lgd_width, height=lgd_height, angle=0, facecolor=color, edgecolor=color)
        except ValueError as error:
            print(f'\nValueError: cannot plot legend for {lbl}: {error}\n')
            exit()
        except Exception as error:
            print(f'\nError: cannot plot legend for {lbl}: {error}\n')
            exit()
            
        ax.add_patch(ellipse)
        ax.text(x+lgd_x_gap, y, lbl, fontsize=graph_params.lgd_font_size)
        row = row + 1

def PlotGrafAnc(args):
    """
    Plot the scatterplot
    """
    
    in_file = args.input_file
    out_file = args.output_file

    # Get graph parameters from GraphParameter class
    graph_params = GraphParameter(args)
    xcol = graph_params.xcol
    ycol = graph_params.ycol
    gwidth = graph_params.width
    gheight = graph_params.height
    xmin = graph_params.xmin
    xmax = graph_params.xmax
    ymin = graph_params.ymin
    ymax = graph_params.ymax

    # Read subject races and get race colors from SubjectPops class
    sbj_pop_obj = SubjectPops(args)
    sbj_pops = sbj_pop_obj.sbj_pops
    pop_colors = sbj_pop_obj.pop_colors
    num_sbjs = sbj_pop_obj.num_sbjs
    num_pops = sbj_pop_obj.num_pops

    fig, ax = plt.subplots(figsize=(gwidth, gheight))
    if args.tight_plot:
        fig, ax = plt.subplots(figsize=(gwidth, gheight), layout='constrained')

    # Read subjects and their GrafAnc scores from the input file
    pf = None
    try:
        pf = pd.read_csv(in_file, sep="\t")
    except Exception as error:
        print(f'\nERROR: Cannot open {in_file}: {error}\n')
        exit()
        
    num_file_sbjs = len(pf[sbj_pop_obj.ga_sbj_col]) # not all will be read
    
    sbjs = []
    
    # Selected populations are plotted to the front
    xvals = []
    yvals = []
    colors = []

    # Other populations are plotted to the background using background color
    xvals0 = []
    yvals0 = []
    colors0 = []

    for i in range(num_file_sbjs):
        sbj = str(pf[sbj_pop_obj.ga_sbj_col][i])
        # Plot only those with Anc IDs from the SubjectPops object
        if sbj not in sbj_pop_obj.sbj_anc_ids.keys():
            continue

        pop = sbj_pop_obj.sbj_pops.get(sbj)
        if pop is None:
            pop = "NOT REPORTED"

        anc = sbj_pop_obj.sbj_anc_ids.get(sbj)
        sup_anc = sbj_pop_obj.sbj_sup_anc_ids.get(sbj)

        if anc is None or sup_anc is None:
            print(f'\nERROR: no Anc ID assigned to subject {sbj}\n')
            exit()

        color = sbj_pop_obj.background_color
        if pop in sbj_pop_obj.pop_colors:
            color = sbj_pop_obj.pop_colors[pop]

        x = pf[xcol][i]
        y = pf[ycol][i]

        if x > xmin and x < xmax and y > ymin and y < ymax:
            sbjs.append(sbj)

            # Those with background color are unselected and to be set to background            
            if color == sbj_pop_obj.background_color:
                colors0.append(color)
                xvals0.append(x)
                yvals0.append(y)
            else:
                colors.append(color)
                xvals.append(x)
                yvals.append(y)

    # Plot the AF line for GD1 vs. GD3
    if not args.hide_EAF:
        e_x, e_y = graph_params.e_x, graph_params.e_y
        f_x, f_y = graph_params.f_x, graph_params.f_y
        a_x, a_y = graph_params.a_x, graph_params.a_y
        line_color = "gray"
        line_width = 0.8
        
        # Plot the EAF triangle for GD1 vs. GD2
        if xcol == "GD1" and ycol == "GD2":
            ax.plot([e_x, a_x, f_x, e_x],[e_y, a_y, f_y, e_y], linewidth=line_width, color=line_color)
                
        # Plot the AF line for GD1 vs. GD3
        if xcol == "GD1" and ycol == "GD3":
            ax.plot([a_x, f_x],[0, 0], linewidth=0.8, color=line_color)

    ShowLegends(ax, graph_params, sbj_pop_obj)

    if not args.nolimits:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

    ax.scatter(xvals0, yvals0, s=graph_params.dot_size, c=colors0, data=pf)
    ax.scatter(xvals, yvals, s=graph_params.dot_size, c=colors, data=pf)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.tick_params(axis='x', labelsize=graph_params.tick_lbl_size)
    ax.tick_params(axis='y', labelsize=graph_params.tick_lbl_size)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%2.1f"))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%2.1f"))
    
    if graph_params.plot_title != "":
        ax.set_title(graph_params.plot_title)

    try:
        plt.savefig(out_file)
    except Exception as error:
        print(f'\nERROR: Cannot save results to {out_file}: {error}\n')
    else:
        print(f'\nSaved scatterplot to {out_file}\n')

def main():
    """
    This script takes the file generated by the C++ program grafanc and plots
    results to a scatterplot, which is saved to the specified output .png file.

    Author: Yumi (Jimmy) Jin
    Version: 0.4.0    
    Date: 04/14/2025
    """
    parser = argparse.ArgumentParser(description="This script plots GrafAnc results to a scatterplot.")
    parser.add_argument("input_file", help="File generated by C++ program grafanc")
    parser.add_argument("output_file", help="The .png file to save the scatterplot")
    parser.add_argument("x_score", help="Score (e.g., GD1) to be plotted on the x-axis")
    parser.add_argument("y_score", help="Score (e.g., GD2) to be plotted on the y-axis")
    parser.add_argument("--sbj_pop_file", help="Tab-delimited file containing self-reported subject races", metavar="")
    parser.add_argument("--sbj_col", help="Header of the subject column", metavar="")
    parser.add_argument("--pop_col", help="Header of the race/ethnicity column", metavar="")
    parser.add_argument("--sub_anc", help="To be colored by subcontinental ancestry group", action="store_true")
    parser.add_argument("--show_ancs", help="Specify Anc IDs to display (comma-delimitted numbers like 100,203)", metavar="")
    parser.add_argument("--pops", help="Populations to be plotted (comma-delimitted numbers like 1,5,3)", metavar="")
    parser.add_argument("--color_file", help="File containing user-specified color list", metavar="")
    parser.add_argument("--bg_color", help="Background color for unselected populations", metavar="")
    parser.add_argument("--color_str", help="List of user-specified colors (comma delimited)", metavar="")
    parser.add_argument("--gw", help="Graph width in inches", metavar="")
    parser.add_argument("--gh", help="Graph height in inches", metavar="")
    parser.add_argument("--xmin", help="Min x value to plot", metavar="")
    parser.add_argument("--xmax", help="Max x value to plot", metavar="")
    parser.add_argument("--ymin", help="Min y value to plot", metavar="")
    parser.add_argument("--ymax", help="Max y value to plot", metavar="")
    parser.add_argument("--lgd_x", help="Legend x position (0=left, 0.5=middle, 1=right)", metavar="")
    parser.add_argument("--lgd_y", help="Legend y position (0=top, 0.5=middle, 1=bottom)", metavar="")
    parser.add_argument("--dot_size", help="Dot size (0.5, 1, 3, etc) of the scatter plot", metavar="")
    parser.add_argument("--title", help="The title of the scatterplot", metavar="")
    parser.add_argument("--tight_plot", help="Generate a tight plot without margins", action="store_true")
    parser.add_argument("--hide_EAF", help="Do not show EAF triangle on GD1-GD2 plot", action="store_true")
    parser.add_argument("--nolimits", help="Let python select axis limits", action="store_true")

    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        print(f"\nERROR: didn't find file {args.input_file}\n")
        exit()

    PlotGrafAnc(args)
    
if __name__ == "__main__":
    main()
    
