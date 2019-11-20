import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None


class MajorRestPeriod(object):
    """
    Calculates the major rest period on a single day (24 hours) of geneactiv accelerometer data.

    """

    def __init__(
        self,
        df,
        fs,
        temperature_threshold=25.0,
        minimum_rest_block=30,
        allowed_rest_break=60,
        minimum_rest_threshold=0.1,
        maximum_rest_threshold=1000.0,
    ):
        """
        Class initialization

        Parameters
        ----------
        df : data frame
            Geneactiv data loaded into a pandas dataframe. Expected duration of 24 hours.
        fs : float
            Sampling rate of the input data.
        temperature_threshold : float
            Lowest temperature for which a data point is considered valid.
        minimum_rest_block : int
            Number of minutes required to consider a rest period valid.
        allowed_rest_break : int
            Number of minutes allowed to interrupt major rest period.
        minimum_rest_threshold : float
            Minimum allowed threshold for determining major rest period.
        maximum_rest_threshold : float
            Maximum allowed threshold for determining major rest period.

        """
        self.df = df
        self.df_angle = None
        self.fs = fs
        self.min_t = temperature_threshold
        self.minimum_rest_block = minimum_rest_block
        self.allowed_rest_break = allowed_rest_break
        self.minimum_rest_threshold = minimum_rest_threshold
        self.maximum_rest_threshold = maximum_rest_threshold
        self.rest_periods = None
        self.mrp = []

    def run(self):
        """
        Runs the algorithm to detect major rest period.

        Returns
        -------
        mrp : list of timestamp
            Major rest period for the given 24 hour period.

        """

        # RELEVANT STREAMS
        df = self.df[["X", "Y", "Z", "T"]]

        # ROLLING 5S MEDIAN
        df = df.rolling(int(5 * self.fs), center=True).median()
        df["z_angle"] = np.arctan(df["Z"] / ((df["X"] ** 2 + df["Y"] ** 2) ** 0.5)) * (
            180.0 / np.pi
        )

        # 5S MEAN
        df = df[["z_angle", "T"]].resample("5s").mean()

        # KEEP ANGLE FOR PLOTTING
        self.df_angle = pd.DataFrame(df["z_angle"]).copy()

        # ABSOLUTE DIFFERENCE
        df["z_angle"] = np.abs(df["z_angle"] - df["z_angle"].shift(1))

        # ROLLING 5MIN MEDIAN
        df_angle = pd.DataFrame(df["z_angle"].rolling(60, center=True).median())
        df_temp = pd.DataFrame(df["T"].rolling(60, center=True).median())

        # CALCULATE THRESHOLD
        thresh = np.min(
            [
                np.max(
                    [
                        np.percentile(df_angle.z_angle.dropna().values, 10) * 15.0,
                        self.minimum_rest_threshold,
                    ]
                ),
                self.maximum_rest_threshold,
            ]
        )

        # APPLY THRESHOLD
        df_angle.z_angle[df_angle.z_angle < thresh] = 0
        df_angle.z_angle[df_angle.z_angle >= thresh] = 1

        # DROP REST PERIODS WITH TEMP TOO LOW
        df_angle.z_angle[df_temp["T"] <= self.min_t] = 1
        df = df_angle.fillna(value=1)

        # DROP REST BLOCKS < MIN
        df["block"] = (df.z_angle.diff().ne(0)).cumsum()
        groups, iter_count = df.groupby(by="block"), 0
        for group in groups:
            iter_count += 1
            if iter_count == 1 or iter_count == len(groups):
                continue
            if (
                group[1]["z_angle"].sum() == 0
                and len(group[1]) < 12 * self.minimum_rest_block
            ):
                df.z_angle[group[1].index[0] : group[1].index[-1]] = 1

        # DROP ACTIVE BLOCKS < MIN
        df["block"] = (df.z_angle.diff().ne(0)).cumsum()
        groups, iter_count = df.groupby(by="block"), 0
        for group in groups:
            iter_count += 1
            if iter_count == 1 or iter_count == len(groups):
                continue
            if (
                len(group[1]) == group[1]["z_angle"].sum()
                and len(group[1]) < 12 * self.allowed_rest_break
            ):
                df.z_angle[group[1].index[0] : group[1].index[-1]] = 0

        # KEEP CANDIDATE REST PERIODS FOR PLOTTING
        self.rest_periods = pd.DataFrame(df.z_angle)

        # LONGEST BLOCK
        df["block"] = (df.z_angle.diff().ne(0)).cumsum()
        best = 1
        for group in df.groupby(by="block"):
            if group[1]["z_angle"].sum() == 0 and len(group[1]) > best:
                best = len(group[1])
                self.mrp = [group[1].index[0], group[1].index[-1]]
        return self.mrp
