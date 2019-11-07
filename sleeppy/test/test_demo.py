import pytest
import time
import os
from importlib import resources
from sleeppy.sleep import SleepPy


@pytest.mark.skip(reason="Run only manually")
def test_demo():
    with resources.path("sleeppy.test.test_data", "demo.bin") as src:
        src = str(src)
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
    test_demo()
