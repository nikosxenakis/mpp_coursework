#!/usr/bin/env python2
# data_analyzer.py
import os
import matplotlib.pyplot as plt


def create_plots(title, x_axis_title, y_axis_title, labels, x_values, y_values, alpha, legent_pos):

    fig, ax = plt.subplots()
    ax.set_xlabel(x_axis_title)
    ax.set_ylabel(y_axis_title)

    colors = ['#35d20a', '#d3390a']

    lines = []

    i = 0
    for y_list in y_values:
        lines.append(plt.plot(x_values, y_list, linewidth=2, color=colors[i]))
        lines.append(plt.plot(x_values, y_list, label=labels[i], linewidth=2, color=colors[i]))
        i = i + 1

    plt.legend(loc=legent_pos)

    fig.savefig('./graphs/' + str(title) + '.eps', format='eps', dpi=1000)

def analyze_data(data_file, image):

    i = 0;
    times = 10
    curr_times = 0
    titles = ()
    processes = ()
    speedup = ()
    sum_run_time = 0
    mean_run_time = 0
    sum_time1 = 0
    mean_time1 = 0
    x = 1

    for line in data_file:
        data = line.split('\t')

        if i == 0:
            titles = data

        if i > 0 and data[0] == str("./resources/" + image + ".pgm"):
            sum_run_time = sum_run_time + float(data[2])
            if int(data[1]) == 1:
                sum_time1 = sum_time1 + float(data[2])

            if x == 1:
                processes = processes + (int(data[1]),)
            if x == 10:
                mean_time1 = sum_time1 / 10
                mean_run_time = sum_run_time / 10
                curr_speedup = mean_time1 / mean_run_time
                speedup = speedup + ( float("%.2f" % curr_speedup ) ,)
                x = 0
                sum_run_time = 0

            x = x + 1
        i = i + 1

    create_plots(
        str(image + "_speedup"),
        titles[1],
        str("speedup"),
        [image, "ideal"],
        processes,
        [speedup, processes],
        0.4,
        "upper left"
    )


def analyze_average_pixel(data_file, path):

    i = 0;

    titles = ()
    iterations = ()
    average = ()

    for line in data_file:
        if i == 0:
            titles = line.split('\t')
        else:
            data = line.split('\t')
            iterations = iterations + (int(data[0]),)
            average = average + ( float("%.3f" % float(data[1])) ,)
        i = i + 1

    create_plots(
        str(path),
        titles[0],
        str(path),
        [titles[1]],
        iterations,
        [average],
        0.4,
        "upper right"
    )


path = './data/results.tsv'

if os.path.exists(path):
    f = open(path, 'r')
    analyze_data(f, "edgenew192x128")
    f = open(path, 'r')
    analyze_data(f, "edgenew256x192")
    f = open(path, 'r')
    analyze_data(f, "edgenew512x384")
    f = open(path, 'r')
    analyze_data(f, "edgenew768x768")

path = './data/'
ext = '.tsv'
filename = 'edgenew192x128_average_pixel'

if os.path.exists(path):
    f = open(path + filename + ext, 'r')
    analyze_average_pixel(f, filename)

filename = 'edgenew256x192_average_pixel'

if os.path.exists(str(path+filename+ext)):
    f = open(str(path+filename+ext), 'r')
    analyze_average_pixel(f, filename)

filename = 'edgenew512x384_average_pixel'

if os.path.exists(str(path+filename+ext)):
    f = open(str(path+filename+ext), 'r')
    analyze_average_pixel(f, filename)

filename = 'edgenew768x768_average_pixel'

if os.path.exists(str(path+filename+ext)):
    f = open(str(path+filename+ext), 'r')
    analyze_average_pixel(f, filename)
