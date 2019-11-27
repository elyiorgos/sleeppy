import pytest
import tempfile
import numpy as np
import pandas as pd
from importlib import resources
from sleeppy.sleep import SleepPy
from sleeppy.utils import make_dummy_data


# INTEGRATION TESTS
def test_sleep_predict():
    # DUMMY DATA
    data, start, end = make_dummy_data()

    # INIT CLASS
    sleep = SleepPy("_", "_", sampling_frequency=20.0)
    sleep.raw_days = [data]

    # PREDICTIONS
    sleep.calculate_sleep_predictions()
    idx_start = sleep.sleep_wake_prediction_days[0].values.argmin()
    idx_end = sleep.sleep_wake_prediction_days[0].values[::-1].argmin()
    test_start = sleep.sleep_wake_prediction_days[0].index[idx_start]
    test_end = sleep.sleep_wake_prediction_days[0].index[::-1][idx_end]
    test_bounds = [test_start.value, test_end.value]

    # SET EXPECTED BOUNDS
    expected_bounds = [i.value for i in [start, end]]

    # DIFFERENCE IN NANOSECONDS
    diffs = np.abs(
        [test_bounds[0] - expected_bounds[0], test_bounds[1] - expected_bounds[1]]
    )

    # MAX DIFFERENCE OF 10min
    max_diff = 600 * 1e9

    # CHECK BOUNDS
    check_bounds = [True if i < max_diff else False for i in diffs]

    # CHECK PREDICTIONS
    awake = sleep.sleep_wake_prediction_days[0].values.sum()
    asleep = len(sleep.sleep_wake_prediction_days[0].values) - awake
    expected_awake = 970
    expected_asleep = 470
    check_predictions = [
        True if awake == expected_awake and asleep == expected_asleep else False
    ]

    assert np.all(check_bounds + check_predictions)


def test_calculate_endpoints():
    # DUMMY DATA
    data, start, end = make_dummy_data()

    # INIT CLASS
    sleep = SleepPy("_", "_", sampling_frequency=20.0)
    sleep.raw_days = [data]
    sleep.major_rest_periods = [[start, end]]

    # PREDICTIONS
    sleep.calculate_sleep_predictions()

    # ENDPOINTS
    sleep.calculate_endpoints()
    test_endpoints = sleep.endpoints[1]
    expected_endpoints = [470, 97.713, 0, 8, 0]
    assert np.array_equal(test_endpoints, expected_endpoints)


@pytest.mark.skip(reason="Run only manually")
def test_visualize():
    # DUMMY DATA
    data, start, end = make_dummy_data()

    # DST DIR
    with resources.path("sleeppy.test", "test_data") as file_path:
        # INIT CLASS WITH 2 DAYS OF DATA
        sleep = SleepPy(
            "dummy/test_report.csv", str(file_path), sampling_frequency=20.0
        )
        sleep.raw_days = [data, data]
        sleep.raw_days_to_plot = [
            data.resample("60s").median(),
            data.resample("60s").median(),
        ]

        # PREDICTIONS
        sleep.calculate_major_rest_periods()
        sleep.calculate_sleep_predictions()
        sleep.calculate_endpoints()

        # RUN VISUALIZATION
        sleep.visualize()


def test_export(output_table):
    # DUMMY DATA
    data, start, end = make_dummy_data()

    # TEMP OUTPUT DIR
    with tempfile.TemporaryDirectory() as out_dir:

        # INIT CLASS WITH 2 DAYS OF DATA
        sleep = SleepPy("dummy/test_report.csv", out_dir, sampling_frequency=20.0)
        sleep.raw_days = [data, data]
        sleep.raw_days_to_plot = [
            data.resample("60s").median(),
            data.resample("60s").median(),
        ]

        # PREDICTIONS
        sleep.calculate_major_rest_periods()
        sleep.calculate_sleep_predictions()
        sleep.calculate_endpoints()

        # EXPORT DATA
        sleep.export()

        # LOAD EXPECTED AND TEST
        expected = pd.read_csv(output_table("expected"))
        test = pd.read_csv(out_dir + "/test_report_output_table.csv")

        assert np.array_equal(expected, test)
