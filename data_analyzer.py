#!/usr/bin/env python2
# data_analyzer.py
import os
import matplotlib.pyplot as plt


def create_plots(title, x_axis_title, y_axis_title, labels, x_values, y_values, alpha, legent_pos):

    fig, ax = plt.subplots()
    ax.set_xlabel(x_axis_title)
    ax.set_ylabel(y_axis_title)

    colors = ['#d3390a', '#35d20a']

    lines = []

    i = 0
    for y_list in y_values:
        lines.append(plt.plot(x_values, y_list, label=labels[i], linewidth=2, color=colors[i]))
        i = i + 1

    plt.legend(loc=legent_pos)

    fig.savefig('./graphs/' + str(title) + '.eps', format='eps', dpi=1000)

def analyze_data(data_file, image):

    i = 0;

    titles = ()
    processes = ()
    speedup = ()
    time1 = 0

    for line in data_file:
        if i == 0:
            titles = line.split('\t')
        else:
            data = line.split('\t')
            if data[0] == str("./resources/" + image):
                if int(data[1]) == 1:
                    time1 = int(data[1])
                processes = processes + (float(data[1]),)
                speedup = speedup + ( float("%.3f" % float(time1/float(data[2]))) ,)
        i = i + 1

    create_plots(
        str(image + "_speedup"),
        titles[0],
        str(image + "_speedup"),
        [titles[1], titles[2]],
        processes,
        [speedup],
        0.4,
        "upper right"
    )


path = './data/results.tsv'

if os.path.exists(path):
    f = open(path, 'r')
    analyze_data(f, "edgenew192x128.pgm")
    f = open(path, 'r')
    analyze_data(f, "edgenew256x192.pgm")
    f = open(path, 'r')
    analyze_data(f, "edgenew512x384.pgm")
    f = open(path, 'r')
    analyze_data(f, "edgenew768x768.pgm")
