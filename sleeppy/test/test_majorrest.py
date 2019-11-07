import pytest
import numpy as np
from sleeppy.utils import make_dummy_data
from sleeppy.majorrest import MajorRestPeriod


def test_mrp():
    # DUMMY DATA
    input_df, start_sleep, end_sleep = make_dummy_data()

    # CALCULATE MRP
    mrp_object = MajorRestPeriod(input_df, 20.0)
    test_mrp = mrp_object.run()
    test_mrp = [i.value for i in test_mrp]

    # SET EXPECTED MRP
    expected_mrp = [i.value for i in [start_sleep, end_sleep]]

    # DIFFERENCE IN NANOSECONDS
    diffs = np.abs([test_mrp[0] - expected_mrp[0], test_mrp[1] - expected_mrp[1]])

    # MAX DIFFERENCE OF 30s
    max_diff = 30 * 1e9

    # CHECK
    check = [True if i < max_diff else False for i in diffs]

    assert np.all(check)
