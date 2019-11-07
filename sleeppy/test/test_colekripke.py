import pytest
from sleeppy.colekripke import ColeKripke
import pandas as pd
import numpy as np


def test_apply_kernel(activity_index_data):
    # DUMMY DATA
    input_ai = pd.read_hdf(activity_index_data, key="test_ai")

    # RUN KERNEL
    ck = ColeKripke(input_ai)
    output_test = ck.apply_kernel()

    # EXPECTED OUTPUT
    output_expected = pd.read_hdf(
        activity_index_data, key="expected_output"
    ).values.flatten()

    assert np.array_equal(output_test, output_expected)


def test_rescore(activity_index_data):
    # DUMMY DATA
    input_array = pd.read_hdf(
        activity_index_data, key="expected_output"
    ).values.flatten()

    # SLEEP WINDOW (10pm to 6am at 20hz plus 0.5hour buffer)
    input_array = input_array[570:1110]

    # TEST OUTPUT
    ck = ColeKripke(input_array)
    ck.predictions = input_array
    output_test = ck.rescore()

    # EXPECTED OUTPUT
    output_expected = pd.read_hdf(
        activity_index_data, key="expected_output_rescored"
    ).values.flatten()
    assert np.array_equal(output_expected, output_test)
