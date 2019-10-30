import time
import os
from sleeppy.sleep import SleepPy


def run_demo():
    src = "test_data/demo.bin"
    dst = input("Please provide a path to a results directory:    ")
    while not os.path.isdir(dst):
        dst = input(
            "\nThe directory you provided does not exist."
            "\nYour input should follow a format similar to /Users/username/Desktop/Results"
            "\nPlease provide a path to a results directory:    "
        )
    st = time.time()

    try:
        SleepPy(input_file=src, results_directory=dst, sampling_frequency=100).run()
    except Exception as e:
        print("Error processing: {}\nError: {}".format(src, e))

    stp = time.time()
    print("total run time: {} minutes".format((stp - st) / 60.0))


if __name__ == "__main__":
    run_demo()
