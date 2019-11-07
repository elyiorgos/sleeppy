import pytest
import numpy as np
import pandas as pd
from sleeppy.utils import (
    activity_index,
    non_overlapping_windows,
    bin2df,
    high_pass_filter,
    load_geneactiv_csv,
    downsample,
    make_dummy_data,
    total_sleep_time,
    percent_time_asleep,
    wake_after_sleep_onset,
    sleep_onset_latency,
    number_of_wake_bouts,
)


def test_activity_index():
    input_array = np.array(
        (
            [[1, 2, 3], [2, 3, 2], [3, 4, 2], [4, 5, 2], [5, 6, 2], [6, 7, 2]],
            [[1, 2, 3], [2, 3, 2], [3, 4, 2], [4, 5, 2], [5, 6, 2], [6, 7, 2]],
            [[1, 2, 3], [2, 3, 2], [3, 4, 2], [4, 5, 2], [5, 6, 2], [6, 7, 2]],
        )
    )
    test_output = activity_index(input_array)
    expected_output = np.array([[1.33333333], [1.33333333], [1.33333333]])
    assert np.allclose(test_output, expected_output)


def test_non_overlapping_windows():
    input_df = pd.DataFrame(np.zeros([100, 3]))
    output = non_overlapping_windows(input_df, 3)
    test_shape = output.shape
    expected_shape = (33, 3, 3)
    assert np.array_equal(test_shape, expected_shape)


def test_bin2df(geneactiv_data):
    input_path = geneactiv_data("bin")
    expected_path = geneactiv_data("csv")
    test_df = bin2df(input_path)
    expected_df = load_geneactiv_csv(expected_path)
    assert np.allclose(test_df, expected_df)


def test_high_pass_filter():
    sub_array = [
        [1, 2, 3],
        [2, 3, 2],
        [3, 4, 2],
        [4, 5, 2],
        [5, 6, 2],
        [6, 7, 2],
        [1, 2, 3],
        [2, 3, 2],
        [3, 4, 2],
        [4, 5, 2],
        [5, 6, 2],
        [6, 7, 2],
        [1, 2, 3],
        [2, 3, 2],
        [3, 4, 2],
        [4, 5, 2],
        [5, 6, 2],
        [6, 7, 2],
        [1, 2, 3],
        [2, 3, 2],
        [3, 4, 2],
        [4, 5, 2],
        [5, 6, 2],
        [6, 7, 2],
        [1, 2, 3],
        [2, 3, 2],
        [3, 4, 2],
        [4, 5, 2],
        [5, 6, 2],
        [6, 7, 2],
        [1, 2, 3],
        [2, 3, 2],
        [3, 4, 2],
        [4, 5, 2],
        [5, 6, 2],
        [6, 7, 2],
    ]
    input_array = np.array((sub_array, sub_array, sub_array))
    sub_array = [
        [-0.14760574, -0.14760574, -0.02994791],
        [0.71247201, 0.71247201, -1.01133389],
        [1.57024711, 1.57024711, -0.99253225],
        [2.42570025, 2.42570025, -0.97354249],
        [3.27881204, 3.27881204, -0.95436413],
        [4.12956311, 4.12956311, -0.93499668],
        [-1.02206592, -1.02206592, 0.08456033],
        [-0.17609447, -0.17609447, -0.89569262],
        [0.66745804, 0.66745804, -0.87575505],
        [1.50857217, 1.50857217, -0.8556265],
        [2.34722848, 2.34722848, -0.83530651],
        [3.18340751, 3.18340751, -0.8147946],
        [-1.98291023, -1.98291023, 0.20590968],
        [-1.15174422, -1.15174422, -0.7731932],
        [-0.32311395, -0.32311395, -0.7521028],
        [0.50296107, 0.50296107, -0.73081865],
        [1.32646133, 1.32646133, -0.70934031],
        [2.14736729, 2.14736729, -0.68766733],
        [-3.03434058, -3.03434058, 0.33420076],
        [-2.21868183, -2.21868183, -0.64373562],
        [-1.40567599, -1.40567599, -0.621476],
        [-0.59534263, -0.59534263, -0.59901996],
        [0.21229869, 0.21229869, -0.57636703],
        [1.0172284, 1.0172284, -0.55351678],
        [-4.18057307, -4.18057307, 0.46953124],
        [-3.38112529, -3.38112529, -0.50722254],
        [-2.58444786, -2.58444786, -0.48377768],
        [-1.79056036, -1.79056036, -0.46013373],
        [-0.99948237, -0.99948237, -0.43629025],
        [-0.21123348, -0.21123348, -0.41224681],
        [-5.42583328, -5.42583328, 0.61199703],
        [-4.64330138, -4.64330138, -0.36355829],
        [-3.86365736, -3.86365736, -0.33891233],
        [-3.08692082, -3.08692082, -0.31406466],
        [-2.31311136, -2.31311136, -0.28901483],
        [-1.54224858, -1.54224858, -0.26376242],
    ]
    expected_output = np.array([sub_array, sub_array, sub_array])
    test_output = high_pass_filter(input_array, 100.0, 0.25, 3)

    assert np.allclose(expected_output, test_output)


def test_downsample():
    # 100hz DUMMY DATA
    df, _, _ = make_dummy_data(freq="10L")

    # DOWNSAMPLE TO 20hz (factor of 5
    df_ds = downsample(df, 5)

    # ENSURE STRUCTURE OF DF REMAINS CONSTANT
    check = [
        np.array_equal(df.columns, df_ds.columns),
        (df.index[0] == df_ds.index[0]),
        (df.index[-1] == df_ds.index[-1]),
    ]
    assert np.all(check)


def test_tst():
    # DUMMY DATA
    predictions = np.empty((1000,), int)
    predictions[::2] = 1.0
    predictions[1::2] = 0.0

    expected_output = 500

    assert expected_output == total_sleep_time(predictions)


def test_waso():
    # DUMMY DATA
    predictions = np.empty((1000,), int)
    predictions[::2] = 1.0
    predictions[1::2] = 0.0

    expected_output = 499

    assert expected_output == wake_after_sleep_onset(predictions)


def test_pta():
    # DUMMY DATA
    predictions = np.empty((1000,), int)
    predictions[::2] = 1.0
    predictions[1::2] = 0.0

    expected_output = 50.00

    assert expected_output == percent_time_asleep(predictions)


def test_nwb():
    # DUMMY DATA
    predictions = np.empty((1000,), int)
    predictions[::2] = 1.0
    predictions[1::2] = 0.0

    expected_output = 499

    assert expected_output == number_of_wake_bouts(predictions)


def test_sol():
    # DUMMY DATA
    predictions = np.empty((1000,), int)
    predictions[::2] = 1.0
    predictions[1::2] = 0.0

    expected_output = 1

    assert expected_output == sleep_onset_latency(predictions)
