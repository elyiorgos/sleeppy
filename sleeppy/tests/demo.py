import time
import os
from sleeppy.sleep import SleepPy
import pandas as pd


def run_demo(binder_demo=False):
    src = __file__.split(".py")[0] + ".bin"
    dst = raw_input("Please provide a path to a results directory:    ")
    while not os.path.isdir(dst):
        if not binder_demo:
            dst = raw_input(
                "\nYour previous entry was not appropriate."
                "\nIt should follow a format similar to /Users/username/Desktop/Results"
                "\nPlease provide a path to a results directory:    "
            )
        else:
            dst = raw_input(
                "\nYour previous entry was not appropriate."
                "\nIt should follow a format similar to /home/jovyan/example_notebook"
                "\nPlease provide a path to a results directory:    "
            )

    st = time.time()
    try:
        if not binder_demo:
            SleepPy(input_file=src, results_directory=dst, sampling_frequency=100, verbose=True)
        else:
            SleepPy(input_file=src, results_directory=dst, sampling_frequency=100, run_config=4, verbose=True)
    except Exception as e:
        print("Error processing: {}\nError: {}".format(src, e))
    stp = time.time()
    print("total run time: {} minutes".format((stp - st) / 60.0))

    print("Checking endpoints...")
    expected = [1, 470, 100, 2, 0, 1]
    obtained = collect_endpoints(dst)
    test_endpoints(expected, obtained)


def collect_endpoints(results_dir):
    src = results_dir + "/demo/sleep_endpoints/sleep_endpoints_summary.csv"
    return pd.read_csv(src).values[0]


def test_endpoints(expected, obtained):
    endpoint_names = [
        "day",
        "total_sleep_time",
        "percent_time_asleep",
        "waso",
        "sleep_onset_latency",
        "number_wake_bouts",
    ]
    errors = 0
    for i in range(len(expected)):
        try:
            assert expected[i] == obtained[i]
        except AssertionError:
            print(
                "Error encountered: endpoint at index {} ({}) "
                "does not match expected output".format(i, endpoint_names[i])
            )
            errors += 1
    if errors > 0:
        print("Total of {} errors encountered".format(errors))
    else:
        print("All tests passed")


if __name__ == "__main__":
    run_demo()
