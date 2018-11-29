#!/usr/bin/env python2
# data_analyzer.py
#B136013
#
import os
import matplotlib.pyplot as plt
import csv
from operator import itemgetter
from matplotlib.ticker import MultipleLocator

def sort_timings():
    path = './data/time_results_timing_test.tsv'

    f = open(path, 'r')
    file_str = ""
    i = 0
    for line in f:
        if i > 0:
            file_str = file_str + str(line)
        i = i + 1
    f.close()
    f = open(path, 'w')
    f.write(file_str)
    f.close()

    reader = csv.reader(open(path), delimiter="\t")
    file_str = ""
    for line in sorted(reader, key=lambda x: int(x[1])):
        i = 0
        for item in line:
            if i > 0:
                file_str = file_str + "\t"
            file_str = file_str + str(item)
            i = i + 1
        file_str = file_str + "\n"

    f = open(path, 'w')
    f.write("Input File\tProcesses Number\tAverage Iteration Time (ms)\n")
    f.write(file_str)
    f.close()


def create_plots(title, x_axis_title, y_axis_title, labels, x_values, y_values, alpha, legent_pos):

    fig, ax = plt.subplots()
    ax.set_xlabel(x_axis_title)
    ax.set_ylabel(y_axis_title)

    colors = ['#003f5c', '#444e86', '#955196', '#dd5182', '#ff6e54', '#ffa600']

    lines = []

    i = 0
    for y_list in y_values:
        x_list = x_values[0:len(y_list)]
        lines.append(plt.plot(x_list, y_list, linewidth=2, color=colors[i]))
        lines.append(plt.plot(x_list, y_list, label=labels[i], linewidth=2, color=colors[i]))
        i = i + 1

    # if title == "speedupBigInput":
    #     ml = MultipleLocator(36)
    #     ax.xaxis.set_minor_locator(ml)
    #     ax.xaxis.grid(which="minor", color='k', linestyle='-.', linewidth=0.7)

    plt.legend(loc=legent_pos)

    fig.savefig('./graphs/' + str(title) + '.eps', format='eps', dpi=1000)


def analyze_average_iteration_interval(data_file, image):
    i = 0;
    curr_times = 0
    processes = ()
    speedup = ()
    run_time = 0
    mean_run_time = 0
    time1 = 0

    for line in data_file:
        data = line.split('\t')

        if i > 0 and data[0] == str("./resources/" + image + ".pgm"):
            if int(data[1]) == 1:
                time1 = float(data[2])
            if int(data[1]) <= 16:
                processes = processes + (int(data[1]),)
                run_time = float(data[2])
                curr_speedup = time1 / run_time
                speedup = speedup + ( float("%.2f" % curr_speedup ) ,)

        i = i + 1

    sorted(speedup)
    return speedup


def analyze_average_iteration(data_file, image):
    i = 0;
    curr_times = 0
    processes = ()
    speedup = ()
    run_time = 0
    mean_run_time = 0
    time1 = 0

    for line in data_file:
        data = line.split('\t')

        if i > 0 and data[0] == str("./resources/" + image + ".pgm"):
            if int(data[1]) == 1:
                time1 = float(data[2])
            if int(data[1]) <= 36:
                processes = processes + (int(data[1]),)
                run_time = float(data[2])
                curr_speedup = time1 / run_time
                speedup = speedup + ( float("%.2f" % curr_speedup ) ,)

        i = i + 1

    sorted(speedup)
    return speedup


def get_big_processes(data_file, image):
    i = 0;
    processes = ()

    for line in data_file:
        data = line.split('\t')

        if i > 0 and data[0] == str("./resources/" + image + ".pgm") and int(data[1]) >= 42:
            processes = processes + (int(data[1]),)
        i = i + 1

    sorted(processes)
    return processes


def analyze_average_big_iteration(data_file, image):
    i = 0;
    speedup = ()
    run_time = 0
    time1 = 0

    for line in data_file:
        data = line.split('\t')

        if i > 0 and data[0] == str("./resources/" + image + ".pgm"):

            if int(data[1]) == 1:
                time1 = float(data[2])

            if int(data[1]) >= 42:
                run_time = float(data[2])
                curr_speedup = time1 / run_time
                speedup = speedup + ( float("%.2f" % curr_speedup ) ,)

        i = i + 1

    sorted(speedup)
    return speedup


def get_time1(data_file, image):
    i = 0;
    curr_times = 0
    titles = ()
    processes = ()
    speedup = ()
    run_time = 0
    mean_run_time = 0
    time1 = 0

    for line in data_file:
        data = line.split('\t')

        if i == 0:
            titles = data

        if i > 0 and data[0] == str("./resources/" + image + ".pgm"):
            if int(data[1]) == 1:
                return float(data[2])

        i = i + 1
    return 0


def get_processes(data_file, image):
    i = 0;
    processes = ()

    for line in data_file:
        data = line.split('\t')

        if i > 0 and data[0] == str("./resources/" + image + ".pgm") and int(data[1]) <= 36:
            processes = processes + (int(data[1]),)
        i = i + 1

    sorted(processes)
    return processes


def create_speedup_plot(processes, speedup1, speedup2, speedup3, speedup4, speedup5):
    create_plots(
        str("speedup"),
        "Number of Processes",
        str("Speedup"),
        ["192x128", "256x192", "512x384", "768x768", "1600x1200", "ideal"],
        processes,
        [speedup1, speedup2, speedup3, speedup4, speedup5, processes],
        0.4,
        "upper left"
    )


def create_speedup_big_input_plot(processes, speedup1, speedup2):
    create_plots(
        str("speedupBigInput"),
        "Number of Processes",
        str("Speedup of big images"),
        ["768x768", "1600x1200", "ideal"],
        processes,
        [speedup1, speedup2, processes],
        0.4,
        "upper left"
    )


def create_speedup_interval_plot(processes, speedup, speedup_interval):
    create_plots(
        str("speedupInterval"),
        "Number of Processes",
        str("Speedup of image 768x768"),
        ["Without Intervals", "With Intervals"],
        processes,
        [speedup, speedup_interval],
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
        str("average_pixel/")+str(path),
        "Iterations",
        str("Average Pixel"),
        [path],
        iterations,
        [average],
        0.4,
        "upper right"
    )

    return average


def parse_average_pixel():
    path = './data/average_pixel_test/'
    ext = '.tsv'
    pixel1 = ()
    pixel2 = ()
    pixel3 = ()
    pixel4 = ()
    pixel5 = ()

    filename = 'edgenew192x128_average_pixel'
    f = open(path + filename + ext, 'r')
    pixel1 = analyze_average_pixel(f, filename)

    filename = 'edgenew256x192_average_pixel'
    f = open(str(path+filename+ext), 'r')
    pixel2 = analyze_average_pixel(f, filename)

    filename = 'edgenew512x384_average_pixel'
    f = open(str(path+filename+ext), 'r')
    pixel3 = analyze_average_pixel(f, filename)

    filename = 'edgenew768x768_average_pixel'
    f = open(str(path+filename+ext), 'r')
    pixel4 = analyze_average_pixel(f, filename)

    filename = 'edgenew1600x1200_average_pixel'
    f = open(str(path+filename+ext), 'r')
    pixel5 = analyze_average_pixel(f, filename)


def parse_running_time():
    path = './data/time_results_timing_test.tsv'

    f = open(path, 'r')
    processes = get_processes(f, "edgenew192x128")

    f = open(path, 'r')
    speedup1 = analyze_average_iteration(f, "edgenew192x128")
    f = open(path, 'r')
    speedup2 = analyze_average_iteration(f, "edgenew256x192")
    f = open(path, 'r')
    speedup3 = analyze_average_iteration(f, "edgenew512x384")
    f = open(path, 'r')
    speedup4 = analyze_average_iteration(f, "edgenew768x768")
    f = open(path, 'r')
    speedup5 = analyze_average_iteration(f, "edgenew1600x1200")

    create_speedup_plot(processes, speedup1, speedup2, speedup3, speedup4, speedup5)


def parse_interval():
    processes = ()
    speedup = ()
    speedup_interval = ()

    path = './data/time_results_timing_test.tsv'
    f = open(path, 'r')
    processes = get_processes(f, "edgenew768x768")
    f = open(path, 'r')
    speedup = analyze_average_iteration_interval(f, "edgenew768x768")

    path = './data/time_results_timing_intervals_test.tsv'
    f = open(path, 'r')
    speedup_interval = analyze_average_iteration_interval(f, "edgenew768x768")

    create_speedup_interval_plot(processes, speedup, speedup_interval)


def parse_big_input():

    path = './data/time_results_timing_test.tsv'

    f = open(path, 'r')
    processes = get_big_processes(f, "edgenew768x768")
    f = open(path, 'r')
    speedup1 = analyze_average_big_iteration(f, "edgenew768x768")
    f = open(path, 'r')
    speedup2 = analyze_average_big_iteration(f, "edgenew1600x1200")

    create_speedup_big_input_plot(processes, speedup1, speedup2)


#sorts the timings in time_results_timing_test.tsv file
sort_timings()

#parses and export the graphs of the average pixel for each image
parse_average_pixel()

#parses and export the graphs of the speedup for each image
parse_running_time()

#parses and export the graphs of the speedup with and without the intervals
parse_interval()

#parses and export the graphs of the speedup for each big image
parse_big_input()
