import time
import os
from sleeppy.sleep import *


def run_demo():
    src = __file__.split('.py')[0] + '.bin'
    st = time.time()
    dst = raw_input("Please provide a path to a results directory:    ")
    while not os.path.isdir(dst):
        dst = raw_input("\nYour previous entry was not appropriate."
                        "\nIt should follow a format similar to /Users/username/Desktop/Results"
                        "\nPlease provide a path to a results directory:    ")
    try:
        SleepPy(input_file=src,
                results_directory=dst,
                sampling_frequency=100)
    except Exception as e:
        print "Error processing: {}\nError: {}".format(src, e)
    stp = time.time()
    print "total run time: {} minutes".format((stp-st)/60.)


if __name__ == "__main__":
    run_demo()
