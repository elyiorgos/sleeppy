import unittest
from sleeppy.colekripke import ColeKripke
import pandas as pd
import numpy as np


class TestColeKripke(unittest.TestCase):
    def test_apply_kernel(self):
        # DUMMY DATA
        input_ai = pd.read_hdf("test_data/test_activity_index.h5", key="test_ai")

        # RUN KERNEL
        ck = ColeKripke(input_ai)
        output_test = ck.apply_kernel()

        # EXPECTED OUTPUT
        output_expected = pd.read_hdf(
            "test_data/test_activity_index.h5", key="expected_output"
        ).values.flatten()

        self.assertTrue(np.array_equal(output_test, output_expected))

    def test_rescore(self):
        # DUMMY DATA
        input_array = pd.read_hdf(
            "test_data/test_activity_index.h5", key="expected_output"
        ).values.flatten()

        # SLEEP WINDOW (10pm to 6am at 20hz plus 0.5hour buffer)
        input_array = input_array[570:1110]

        # TEST OUTPUT
        ck = ColeKripke(input_array)
        ck.predictions = input_array
        output_test = ck.rescore()

        # EXPECTED OUTPUT
        output_expected = pd.read_hdf(
            "test_data/test_activity_index.h5", key="expected_output_rescored"
        ).values.flatten()
        self.assertTrue(np.array_equal(output_expected, output_test))


if __name__ == "__main__":
    unittest.main()
