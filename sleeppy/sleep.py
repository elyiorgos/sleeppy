import pandas as pd
import os
import numpy as np
import matplotlib

try:
    matplotlib.use("TkAgg")
except Exception as e:
    print "Error: could not use TkAgg as backend"
    pass
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import struct
from bitstring import BitArray
from scipy import signal

sns.set()
pd.options.mode.chained_assignment = None
__all__ = ["SleepPy", "ColeKripke", "band_pass_filter", "activity_index", "bin2df"]


class SleepPy:
    """
    Processes raw GeneActiv accelerometer data from the wrist to determine various sleep metrics and endpoints.

    *** Support for data from other wrist worn accelerometers will be added in the future. To maximize generality an
    input format will be specified. As long as the input specifications are met, the package should produce the proper
    results. ***

    The sleep window detection and wear detection functions of this package are based off of the following papers, as
    well as their implementation in the R package GGIR:

    van Hees V, Fang Z, Zhao J, Heywood J, Mirkes E, Sabia S, Migueles J (2019). GGIR: Raw Accelerometer Data Analysis.
    doi: 10.5281/zenodo.1051064, R package version 1.9-1, https://CRAN.R-project.org/package=GGIR.

    van Hees V, Fang Z, Langford J, Assah F, Mohammad Mirkes A, da Silva I, Trenell M, White T, Wareham N, Brage S (2014).
    'Autocalibration of accelerometer data or free-living physical activity assessment using local gravity and
    temperature: an evaluation on four continents.' Journal of Applied Physiology, 117(7), 738-744.
    doi: 10.1152/japplphysiol.00421.2014, https://www.physiology.org/doi/10.1152/japplphysiol.00421.2014

    van Hees V, Sabia S, Anderson K, Denton S, Oliver J, Catt M, Abell J, Kivimaki M, Trenell M, Singh-Maoux A (2015).
    'A Novel, Open Access Method to Assess Sleep Duration Using a Wrist-Worn Accelerometer.' PloS One, 10(11).
    doi: 10.1371/journal.pone.0142533, http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0142533.

    The detection of sleep and wake states uses a heuristic model based on the algorithm described in:

    Cole, R.J., Kripke, D.F., Gruen, W.'., Mullaney, D.J., & Gillin, J.C. (1992). Automatic sleep/wake identification
    from wrist activity. Sleep, 15 5, 461-9.

    The activity index feature is based on the index described in:

    Bai J, Di C, Xiao L, Evenson KR, LaCroix AZ, Crainiceanu CM, et al. (2016) An Activity Index for Raw Accelerometry
    Data and Its Comparison with Other Activity Metrics. PLoS ONE 11(8): e0160644.
    https://doi.org/10.1371/journal.pone.0160644

    The test data provided for this package is available from:

    Vincent van Hees, Sarah Charman, & Kirstie Anderson. (2018). Newcastle polysomnography and accelerometer data
    (Version 1.0) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.1160410



    """

    def __init__(
        self,
        input_file,
        results_directory,
        sampling_frequency,
        start_buffer="0s",
        stop_buffer="0s",
        start_time="",
        stop_time="",
    ):
        """
        Class initialization.

        :param input_file: full path to the file to be processed
        :param results_directory: full path to the directory where the results should be saved
        :param sampling_frequency: sampling frequency of the GeneActiv data to be processed
        :param start_buffer: number of seconds to ignore from the beginning of the recording
        :param stop_buffer: number of seconds to ignore from the end of the recording
        :param start_time: time_stamp from which to start looking at data (str) format: "%Y-%m-%d %H:%M:%S:%f"
        "param stop time: time stamp at which to stop looking at data (str) format: "%Y-%m-%d %H:%M:%S:%f"
        """
        self.src = input_file  # save input location
        self.dst = results_directory  # save output location
        self.src_name = input_file.split("/")[-1][0:-4]  # save naming convention
        self.sub_dst = (
            results_directory + "/" + self.src_name
        )  # create output directory
        self.fs = sampling_frequency  # save sampling frequency
        self.window_size = 60  # define window size in seconds
        self.band_pass_cutoff = (
            0.25,
            12.0,
        )  # define the cutoffs for the band pass filter
        self.major_rest_periods = []  # initialize a list to save the major rest periods
        self.start_buffer = start_buffer
        self.stop_buffer = stop_buffer
        self.start_time = start_time
        self.stop_time = stop_time
        self.run()  # run the package

    def run(self):
        """
        Runs the full package on the provided file.

        """
        try:
            os.mkdir(self.sub_dst)  # set up output directory
        except OSError:
            pass
        # split the data into 24 hour periods
        if ".bin" in self.src:
            self.split_days_geneactiv_bin()
        elif ".csv" in self.src:
            self.split_days_geneactiv_csv()

        # extract the activity index feature
        self.extract_activity_index()

        # run wear/on-body detection
        self.wear_detection()

        # run major rest period detection
        self.major_rest_period()

        # run sleep wake predictions on the major rest period
        self.sleep_wake_predict()

        # calculate endpoints based on the above predictions
        self.calculate_endpoints()

        # generates visual reports
        self.visualize_results()

    def split_days_geneactiv_csv(self):
        """
        Splits the GeneActiv accelerometer data into 24 hour chunks, defined from noon to noon.

        """
        try:
            os.mkdir(self.sub_dst + "/raw_days")  # set up output directory
        except OSError:
            pass
        # load data and fix time_stamps
        data = pd.read_csv(
            self.src,
            index_col=0,
            skiprows=100,
            header=None,
            names=["Time", "X", "Y", "Z", "LUX", "Button", "T"],
            usecols=["Time", "X", "Y", "Z", "T"],
            dtype={
                "Time": object,
                "X": np.float64,
                "Y": np.float64,
                "Z": np.float64,
                "LUX": np.int64,
                "Button": bool,
                "T": np.float64,
            },
            low_memory=False,
        )
        data.index = pd.to_datetime(data.index, format="%Y-%m-%d %H:%M:%S:%f").values

        # remove any specified time periods from the beginning and end of the file
        data = data.loc[
            data.index[0]
            + pd.Timedelta(self.start_buffer) : data.index[-1]
            - pd.Timedelta(self.stop_buffer)
        ]

        # cut to defined start and end times if specified
        if self.start_time and self.stop_time:
            self.start_time = pd.to_datetime(
                self.start_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            self.stop_time = pd.to_datetime(
                self.stop_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            data = data.loc[self.start_time : self.stop_time]
        elif self.start_time:
            self.start_time = pd.to_datetime(
                self.start_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            data = data.loc[self.start_time :]
        elif self.stop_time:
            self.stop_time = pd.to_datetime(
                self.stop_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            data = data.loc[: self.stop_time]

        # split data into days from noon to noon
        days = data.groupby(pd.Grouper(level=0, freq="24h", base=12))

        # iterate through days keeping track of the day
        count = 0
        for day in days:
            # save each 24 hour day separately if there's enough data to analyze
            df = day[1].copy()
            available_hours = (len(df) / float(self.fs)) / 3600.0
            if available_hours >= 6:
                count += 1
                dst = "/raw_days/{}_day_{}.h5".format(
                    self.src_name, str(count).zfill(2)
                )
                df.to_hdf(self.sub_dst + dst, key="raw_geneactiv_data_24hr", mode="w")
        return

    def split_days_geneactiv_bin(self):
        """
        Splits the GeneActiv accelerometer data into 24 hour chunks, defined from noon to noon.

        """
        try:
            os.mkdir(self.sub_dst + "/raw_days")  # set up output directory
        except OSError:
            pass
        # load data and fix time_stamps
        data = bin2df(self.src)

        # remove any specified time periods from the beginning and end of the file
        data = data.loc[
            data.index[0]
            + pd.Timedelta(self.start_buffer) : data.index[-1]
            - pd.Timedelta(self.stop_buffer)
        ]

        # cut to defined start and end times if specified
        if self.start_time and self.stop_time:
            self.start_time = pd.to_datetime(
                self.start_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            self.stop_time = pd.to_datetime(
                self.stop_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            data = data.loc[self.start_time : self.stop_time]
        elif self.start_time:
            self.start_time = pd.to_datetime(
                self.start_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            data = data.loc[self.start_time :]
        elif self.stop_time:
            self.stop_time = pd.to_datetime(
                self.stop_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            data = data.loc[: self.stop_time]

        # split data into days from noon to noon
        days = data.groupby(pd.Grouper(level=0, freq="24h", base=12))

        # iterate through days keeping track of the day
        count = 0
        for day in days:
            # save each 24 hour day separately if there's enough data to analyze
            df = day[1].copy()
            available_hours = (len(df) / float(self.fs)) / 3600.0
            if available_hours >= 6:
                count += 1
                dst = "/raw_days/{}_day_{}.h5".format(
                    self.src_name, str(count).zfill(2)
                )
                df.to_hdf(self.sub_dst + dst, key="raw_geneactiv_data_24hr", mode="w")
        return

    def extract_activity_index(self):
        """
        Calculates the activity index feature on each 24 hour day.

        """
        try:
            os.mkdir(self.sub_dst + "/activity_index_days")  # set up output directory
        except OSError:
            pass
        count = 0

        # get days
        days = sorted(
            [
                self.sub_dst + "/raw_days/" + i
                for i in os.listdir(self.sub_dst + "/raw_days/")
            ]
        )
        for day in days:
            count += 1

            # load data
            df = pd.read_hdf(day)
            activity = []
            header = ["Time", "activity_index"]
            idx = 0
            window = int(self.window_size * self.fs)
            incrementer = int(self.window_size * self.fs)

            # iterate through windows
            while idx < len(df) - incrementer:
                # preprocessing: BP Filter
                temp = df[["X", "Y", "Z"]].iloc[idx : idx + window]
                start_time = temp.index[0]
                temp.index = range(len(temp.index))  # reset index
                temp = band_pass_filter(
                    temp, self.fs, bp_cutoff=self.band_pass_cutoff, order=3
                )

                # activity index extraction
                bp_channels = [
                    i for i in temp.columns.values[1:] if "bp" in i
                ]  # band pass filtered channels
                activity.append(
                    [
                        start_time,
                        activity_index(temp, channels=bp_channels).values[0][0],
                    ]
                )
                idx += incrementer

            # save data
            activity = pd.DataFrame(activity)
            activity.columns = header
            activity.set_index("Time", inplace=True)
            dst = "/activity_index_days/{}_activity_index_day_{}.h5".format(
                self.src_name, str(count).zfill(2)
            )
            activity.to_hdf(
                self.sub_dst + dst, key="activity_index_data_24hr", mode="w"
            )

    def wear_detection(self):
        """
        Runs wear detection per each 24 hour chunk, saves both the raw predictions of wear, as well as the
        rescored predictions.

        """
        try:
            os.mkdir(self.sub_dst + "/wear_detection")  # set up output directory
        except OSError:
            pass
        count = 0

        # get days
        days = sorted(
            [
                self.sub_dst + "/raw_days/" + i
                for i in os.listdir(self.sub_dst + "/raw_days/")
            ]
        )
        for day in days:
            df = pd.read_hdf(day)[["X", "Y", "Z"]]
            count += 1

            # get std based classification criteria
            df_std = self.roll_std_60_minute(df)
            df_std[df_std >= 0.013] = 1
            df_std[df_std < 0.013] = 0
            df_std = df_std.sum(axis=1)

            # get range based classification criteria
            df_range = self.roll_max_range_60_minute(df)
            df_range[df_range >= 0.15] = 1
            df_range[df_range < 0.15] = 0
            df_range = df_range.sum(axis=1)

            # classify
            df_wear = pd.DataFrame(df_std.copy()) * 0 + 1
            df_wear.columns = ["wear"]
            for i in range(len(df_wear)):
                if df_range.ix[i] <= 1 or df_std.ix[i] <= 1:
                    df_wear.ix[i] = 0

            # save before rescoring
            df_wear.to_hdf(
                self.sub_dst
                + "/wear_detection/wear_detection_day_{}.hdf".format(
                    str(count).zfill(2)
                ),
                key="wear_detection_24hr",
                mode="w",
            )

            # apply rescoring
            df_wear = self.rescore(df_wear)
            df_wear = self.rescore(df_wear)
            df_wear = self.rescore(df_wear)
            if count == len(days):
                df_wear = self.rescore_last_day(df_wear)

            # save post rescoring
            df_wear.to_hdf(
                self.sub_dst
                + "/wear_detection/wear_detection_rescored_day_{}.hdf".format(
                    str(count).zfill(2)
                ),
                key="wear_detection_rescored_24hr",
                mode="w",
            )

    def major_rest_period(self):
        """
        Determines the major rest period,
        """
        try:
            os.mkdir(self.sub_dst + "/major_rest_period")  # set up output directory
        except OSError:
            pass
        count = 0
        mrps = []
        header = ["day", "major_rest_period", "available_hours"]

        # get days
        days = sorted(
            [
                self.sub_dst + "/raw_days/" + i
                for i in os.listdir(self.sub_dst + "/raw_days/")
            ]
        )
        for day in days:
            df = pd.read_hdf(day)[["X", "Y", "Z"]]
            available_hours = (len(df) / float(self.fs)) / 3600.0
            count += 1

            # process data
            df = df.rolling(int(5 * self.fs)).median()  # run rolling median 5 second
            df = np.arctan(df["Z"] / ((df["X"] ** 2 + df["Y"] ** 2) ** 0.5)) * (
                180.0 / np.pi
            )  # get angle
            df = df.resample("5s").mean().fillna(0)  # get 5 second average

            # save intermediate data for plotting
            df.to_hdf(
                self.sub_dst
                + "/major_rest_period/5_second_average_arm_angle_day_{}.hdf".format(
                    str(count).zfill(2)
                ),
                key="arm_angle_data_24hr",
                mode="w",
            )

            df = np.abs(df - df.shift(1))  # get absolute difference
            df = self.roll_med(df, 60)  # run rolling median 5 minute

            thresh = (
                np.percentile(df.Data.dropna().values, 10) * 15.0
            )  # calculate threshold
            df[df < thresh] = 0  # apply threshold
            df[df >= thresh] = 1  # apply threshold

            # drop rest blocks < 30 minutes
            df["block"] = (df.Data.diff().ne(0)).cumsum()
            for group in df.groupby(by="block"):
                if group[1]["Data"].sum() == 0 and len(group[1]) < 360:
                    df.Data[group[1].index[0] : group[1].index[-1]] = 1

            # drop active blocks < 60 minutes
            df["block"] = (df.Data.diff().ne(0)).cumsum()
            for group in df.groupby(by="block"):
                if len(group[1]) == group[1]["Data"].sum() and len(group[1]) < 720:
                    df.Data[group[1].index[0] : group[1].index[-1]] = 0

            # get longest block
            df["block"] = (df.Data.diff().ne(0)).cumsum()
            best = 0
            mrp = []
            for group in df.groupby(by="block"):
                if group[1]["Data"].sum() == 0 and len(group[1]) > best:
                    best = len(group[1])
                    mrp = [group[1].index[0], group[1].index[-1] + pd.Timedelta("5m")]

            # save predictions
            df.drop(columns=["block"], inplace=True)
            df.to_hdf(
                self.sub_dst
                + "/major_rest_period/rest_periods_day_{}.hdf".format(
                    str(count).zfill(2)
                ),
                key="rest_period_data_24hr",
                mode="w",
            )

            mrps.append([count, mrp, available_hours])

        # aggregate and save the major rest period for each day
        mrps = pd.DataFrame(mrps)
        mrps.columns = header
        mrps.set_index("day", inplace=True)
        dst = "/major_rest_period/{}_major_rest_periods.csv".format(self.src_name)
        mrps.to_csv(self.sub_dst + dst)

    def sleep_wake_predict(self):
        """
        Run sleep wake prediction based on the activity index feature.
        """
        try:
            os.mkdir(
                self.sub_dst + "/sleep_wake_predictions"
            )  # set up output directory
        except OSError:
            pass
        count = 0

        # get days
        days = sorted(
            [
                self.sub_dst + "/activity_index_days/" + i
                for i in os.listdir(self.sub_dst + "/activity_index_days/")
            ]
        )
        for day in days:
            count += 1
            df = pd.read_hdf(day)

            # run the sleep wake predictions
            ck = ColeKripke(df.activity_index)
            df["sleep_predictions"] = ck.predict()

            # save predictions
            df.drop(inplace=True, columns=["activity_index"])
            df.to_hdf(
                self.sub_dst
                + "/sleep_wake_predictions/sleep_wake_day_{}.hdf".format(
                    str(count).zfill(2)
                ),
                key="sleep_wake_data_24hr",
                mode="w",
            )

    def calculate_endpoints(self):
        """
        Calculate the metrics and endpoints of interest.
        """
        try:
            os.mkdir(self.sub_dst + "/sleep_endpoints")  # set up output directory
        except OSError:
            pass
        count = 0

        # get days
        days = sorted(
            [
                self.sub_dst + "/sleep_wake_predictions/" + i
                for i in os.listdir(self.sub_dst + "/sleep_wake_predictions/")
            ]
        )

        # get major rest periods for each day
        mrps = pd.read_csv(
            self.sub_dst
            + "/major_rest_period/{}_major_rest_periods.csv".format(self.src_name),
            parse_dates=True,
            index_col="day",
        )
        endpoints = []
        for day in days:
            count += 1
            df = pd.read_hdf(day)
            # get and format times
            times = mrps.loc[count].major_rest_period
            try:
                idt = times.index("[T")
                times = times[: idt + 1] + "pd." + times[idt + 1 :]
                idt = times.index(", ")
                times = times[: idt + 2] + "pd." + times[idt + 2 :]
                times = eval(times)
                df = df.loc[times[0] : times[1]]
            except ValueError:
                pass

            # get total sleep time
            tst = len(df) - sum(df.values)

            # get percent time asleep
            pct_time_sleep = 100.0 * (len(df) - sum(df.values)) / float(len(df))

            # get wake after sleep onset
            waso = df.loc[df.idxmin()[0] :]
            waso = waso.sum()[0]

            # get sleep onset latency
            sleep_onset_lat = (df.idxmin()[0] - df.index[0]).total_seconds() / 60.0

            # number of wake bouts
            num_wake_bouts = 0
            wake_bout_df = df.copy()
            wake_bout_df["block"] = (
                wake_bout_df.sleep_predictions.diff().ne(0)
            ).cumsum()
            for group in wake_bout_df.groupby(by="block"):
                if group[1]["sleep_predictions"].sum() > 0:
                    num_wake_bouts += 1
            endpoints.append(
                [
                    int(count),
                    int(tst[0]),
                    int(np.round(pct_time_sleep[0])),
                    int(waso),
                    int(sleep_onset_lat),
                    int(num_wake_bouts),
                ]
            )

        # build and save output dataframe
        hdr = [
            "day",
            "total_sleep_time",
            "percent_time_asleep",
            "waso",
            "sleep_onset_latency",
            "number_wake_bouts",
        ]
        endpoints = pd.DataFrame(endpoints)
        endpoints.columns = hdr
        endpoints.set_index(endpoints.day, inplace=True)
        endpoints.drop(columns="day", inplace=True)
        endpoints.to_csv(self.sub_dst + "/sleep_endpoints/sleep_endpoints_summary.csv")

    def roll_med(self, df, num_samples):
        """
        Calculate the rolling median of a pandas dataframe.

        :param df: pandas dataframe
        :param num_samples: number of samples to include in the window
        :return: pandas dataframe containing the rolling median values.
        """

        # initialize indexer and rolling median list
        idx = 0
        med = []
        while idx < len(df) - num_samples:
            med.append(
                [df.index[idx], df.iloc[idx : idx + num_samples].median()]
            )  # get start index, std value
            idx += 1
        # format data frame
        df = pd.DataFrame(med, columns=["Time", "Data"])
        df.set_index(df.Time, inplace=True)
        df.drop(inplace=True, columns="Time")
        return df

    def roll_std_60_minute(self, df):
        """
        Calculate the rolling 60 minute standard deviation of a pandas dataframe.

        :param df: pandas dataframe
        :return: pandas dataframe containing the rolling std values
        """
        # initialize indexer and rolling std list
        idx = 0
        rstd = []

        # calculate std for all windows
        while idx < len(df) - int(900 * self.fs):  # run until we reach the end
            xyz = (
                df.ix[idx : idx + int(3600 * self.fs)].std().values
            )  # save std in x y and z
            rstd.append(
                [df.index[idx], xyz[0], xyz[1], xyz[2]]
            )  # get start index of window, std values
            idx += int(900 * self.fs)  # increment indexer by 15 minutes

        # format dataframe
        df = pd.DataFrame(rstd, columns=["Time", "X", "Y", "Z"])
        df.set_index(df.Time, inplace=True)
        df.drop(inplace=True, columns="Time")
        return df

    def roll_max_range_60_minute(self, df):
        """
        Calculate the rolling 60 minute max range of a pandas dataframe.

        :param df: pandas dataframe
        :return: pandas dataframe containing the rolling std values.
        """
        # initialize indexer and rolling range list
        idx = 0
        rr = []

        # calculate range for all windows
        while idx < len(df) - int(900 * self.fs):  # run until we reach the end
            xyz = (
                df.ix[idx : idx + int(3600 * self.fs)].max().values
                - df.ix[idx : idx + int(3600 * self.fs)].min().values
            )  # save range in x y z
            rr.append(
                [df.index[idx], xyz[0], xyz[1], xyz[2]]
            )  # get start index of window, range values
            idx += int(900 * self.fs)  # increment indexer by 15 minutes

        # format dataframe
        df = pd.DataFrame(rr, columns=["Time", "X", "Y", "Z"])
        df.set_index(df.Time, inplace=True)
        df.drop(inplace=True, columns="Time")
        return df

    def visualize_results(self):
        """
        Generates reports to visualize endpoint summary and day to day endpoint behaviors.
        """
        try:
            os.mkdir(self.sub_dst + "/reports")  # set up output directory
        except OSError:
            pass
        # raw (all 3 axes or vmag-1?)
        rdays = sorted(
            [
                self.sub_dst + "/raw_days/" + i
                for i in os.listdir(self.sub_dst + "/raw_days/")
            ]
        )
        # wear (no rescoring)
        wdays = sorted(
            [
                self.sub_dst + "/wear_detection/" + i
                for i in os.listdir(self.sub_dst + "/wear_detection/")
                if "rescored" not in i and ".hdf" in i
            ]
        )
        # wear (with rescoring)
        wdays_re = sorted(
            [
                self.sub_dst + "/wear_detection/" + i
                for i in os.listdir(self.sub_dst + "/wear_detection/")
                if "rescored" in i
            ]
        )
        # major rest (arm angle)
        mrdays_aa = sorted(
            [
                self.sub_dst + "/major_rest_period/" + i
                for i in os.listdir(self.sub_dst + "/major_rest_period/")
                if "angle" in i
            ]
        )
        # major rest (periods)
        mrdays_rp = sorted(
            [
                self.sub_dst + "/major_rest_period/" + i
                for i in os.listdir(self.sub_dst + "/major_rest_period/")
                if "angle" not in i and "hdf" in i
            ]
        )
        # activity index (full 24 hours)
        aidays = sorted(
            [
                self.sub_dst + "/activity_index_days/" + i
                for i in os.listdir(self.sub_dst + "/activity_index_days/")
            ]
        )
        # sleep wake (full 24 hours)
        swdays = sorted(
            [
                self.sub_dst + "/sleep_wake_predictions/" + i
                for i in os.listdir(self.sub_dst + "/sleep_wake_predictions/")
            ]
        )
        # endpoints (graphs/charts per day)
        endpoints = pd.read_csv(
            self.sub_dst + "/sleep_endpoints/sleep_endpoints_summary.csv",
            index_col="day",
        )

        days = range(0, len(rdays))
        for day in days:
            # read the raw data, downsample for plotting
            raw = pd.read_hdf(rdays[day])  # 10 millisecond period
            raw = raw.resample("60s").median()

            # get shared index
            idx = pd.date_range(
                start=raw.index[0].replace(hour=12, minute=0, second=0, microsecond=0),
                periods=1440,
                freq="60s",
            )
            raw = raw.reindex(idx, fill_value=float("nan"))

            # read the wear data, resample and match index with the raw data
            wear = pd.read_hdf(wdays[day])  # 15 minute period
            wear[wear == 0] = float("nan")
            wear = wear.resample("60s").ffill()
            wear = wear.reindex(idx, fill_value=float("nan"))

            # read the wear data with rescoring, resample and match the raw index
            wear_re = pd.read_hdf(wdays_re[day])  # 15 minute period
            wear_re[wear_re == 0] = float("nan")
            wear_re = wear_re.resample("60s").ffill()
            wear_re = wear_re.reindex(idx, fill_value=float("nan"))

            # read the arm angle data, resample and match the raw index
            angle = pd.read_hdf(mrdays_aa[day])  # 5 second period
            angle = angle.resample("60s").max()
            angle = angle.reindex(idx, fill_value=float("nan"))

            # read the major rest period data, resample and match the raw index
            periods = pd.read_hdf(mrdays_rp[day])  # 5 second period
            periods[periods == 1] = float("nan")
            periods[periods == 0] = 1
            periods = periods.resample("60s").max()
            periods = periods.reindex(idx, fill_value=float("nan"))

            # read the acvitity index data, resample and match the raw index
            aindex = pd.read_hdf(aidays[day])  # 1 minute period
            aindex = aindex.resample("60s").max()
            aindex = aindex.reindex(idx, fill_value=float("nan"))

            # read the sleep wake predictions, resample and match the raw index
            swake = pd.read_hdf(swdays[day])  # 1 minute period
            swake[swake == 0] = float("nan")
            swake = swake.resample("60s").max()
            swake = swake.reindex(idx, fill_value=float("nan"))

            # build a dataframe for plotting certain data streams as straight lines
            df = swake.copy()
            df.columns = ["wake"]
            df["rest periods"] = periods.values - 0.05
            df["on body"] = wear.values - 0.1
            df["on body(rescore)"] = wear_re.values - 0.15
            swake, wear, wear_re, periods = [], [], [], []

            # get day endpoints for plotting of table
            t_labels = (
                "Total Sleep Time(minutes)",
                "Percent Time Asleep",
                "Wake After Sleep Onset(minutes)",
                "Sleep Onset Latency(minutes)",
                "Number of Wake Bouts",
            )
            t_vals = [endpoints.loc[day + 1].values]

            # plotting
            fig, (axt, ax0, ax1, ax2, ax3) = plt.subplots(5, 1, figsize=(30, 10))
            plt.suptitle(
                "Visual Report for Source: {}\nDay: {}\nDate: {}".format(
                    self.src_name, day + 1, idx[0].date()
                ),
                fontsize=25,
            )
            hours = mdates.HourLocator(interval=1)
            h_fmt = mdates.DateFormatter("%H:%M")
            all_axes = (ax0, ax1, ax2, ax3)

            # plot table
            tbl = axt.table(
                cellText=t_vals,
                colLabels=t_labels,
                cellLoc="center",
                rowLoc="center",
                loc="center",
                fontsize=20,
            )
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(24)
            tbl.scale(1.1, 2.4)
            axt.axis("off")

            # plot raw
            raw[["X", "Y", "Z"]].plot(ax=ax0, lw=1).legend(
                bbox_to_anchor=(0, 1), fontsize=20
            )
            ax0.set_ylabel("")
            ax0.set_xlabel("")

            # plot activity index
            aindex.plot(ax=ax1, lw=1, color="#6fc276").legend(
                labels=["activity"], bbox_to_anchor=(0, 0.75), fontsize=20
            )
            ax1.set_ylabel("")
            ax1.set_xlabel("")

            # plot arm angle
            angle.plot(ax=ax2, lw=1, color="#b36ff6").legend(
                labels=["arm angle"], bbox_to_anchor=(0, 0.75), fontsize=20
            )
            ax2.set_ylabel("")
            ax2.set_xlabel("")

            # plot dataframe of 4 streams
            df.plot(ax=ax3, lw=8, x_compat=True).legend(
                bbox_to_anchor=(0, 1.3), fontsize=20
            )
            ax3.set_ylabel("")
            ax3.set_xlabel("")

            # plot formatting
            plt.draw()
            count = 0
            for ax in all_axes:
                count += 1
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.spines["bottom"].set_visible(False)
                ax.spines["left"].set_visible(False)
                ax.grid(False)
                if count < 4:
                    ax.get_xaxis().set_ticks([])
                ax.get_yaxis().set_ticks([])
            ax3.xaxis.set_major_locator(hours)
            ax3.xaxis.set_major_formatter(h_fmt)
            plt.subplots_adjust(wspace=0, hspace=0)
            fig.autofmt_xdate()
            for tick in ax3.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
            plt.savefig(
                self.sub_dst + "/reports/Visual_Results_Day_{}.pdf".format(day + 1)
            )
            plt.close()

        # generate a summary plot from endpoint data
        fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, 1, figsize=(12, 12))
        plt.suptitle("Summary Report for Source: {}".format(self.src_name), fontsize=16)
        all_axes = (ax0, ax1, ax2, ax3, ax4)
        ylabels = [
            "Total Sleep\nTime(min)\nMean: {}".format(
                int(np.round(endpoints.total_sleep_time.mean()))
            ),
            "Percent Time\nAsleep\nMean: {}".format(
                int(np.round(endpoints.percent_time_asleep.mean()))
            ),
            "Wake After\nSleep Onset(min)\nMean: {}".format(
                int(np.round(endpoints.waso.mean()))
            ),
            "Sleep Onset\nLatency(min)\nMean: {}".format(
                int(np.round(endpoints.sleep_onset_latency.mean()))
            ),
            "Number of\nWake Bouts\nMean: {}".format(
                int(np.round(endpoints.number_wake_bouts.mean()))
            ),
        ]

        # plot total sleep time
        endpoints.total_sleep_time.plot.bar(ax=ax0, title="")
        # plot percent time asleep
        endpoints.percent_time_asleep.plot.bar(ax=ax1, title="")
        # plot wake after sleep onset
        endpoints.waso.plot.bar(ax=ax2, title="")
        # plot sleep onset latency
        endpoints.sleep_onset_latency.plot.bar(ax=ax3, title="")
        # plot the number of wake bouts
        endpoints.number_wake_bouts.plot.bar(ax=ax4, title="")
        # plot formatting
        count = 0
        for ax in all_axes:
            count += 1
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.grid(False)
            ax.set_ylabel(ylabels[count - 1], rotation=0, fontsize=12, labelpad=50)
            if count < 5:
                ax.set_xlabel("")
                ax.get_xaxis().set_ticks([])
                ax.get_yaxis().set_ticks([])
            else:
                ax.set_xlabel("Day", fontsize=20)
                ax.get_yaxis().set_ticks([])
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
            for p in ax.patches:
                ax.annotate(
                    np.round(p.get_height(), decimals=2),
                    (p.get_x() + p.get_width() / 2.0, 0),
                    ha="center",
                    va="center",
                    xytext=(0, 10),
                    textcoords="offset points",
                    fontweight="bold",
                )
        plt.subplots_adjust(wspace=0, hspace=0.01)
        plt.xticks(fontsize=20)
        plt.draw()
        plt.savefig(self.sub_dst + "/reports/Summary_Report.pdf")
        plt.close()

    def rescore(self, df):
        """
        Rescores wear detection data saved in a pandas dataframe.

        :param df: pandas dataframe containing wear/nonwear predictions
        :return: rescored pandas dataframe of wear/nonwear predictions
        """
        # group classifications into wear and nonwear blocks
        df["block"] = (df.wear.diff().ne(0)).cumsum()
        blocks = list(df.groupby("block"))

        # iterate through blocks
        for i in range(1, len(blocks) - 1):
            wear = blocks[i][1]["wear"].values[
                0
            ]  # get whether or not the block is wear
            if wear:
                # get hour lengths of the previous, current, and next blocks
                prev, current, post = (
                    len(blocks[i - 1][1]) * 0.25,
                    len(blocks[i][1]) * 0.25,
                    len(blocks[i + 1][1]) * 0.25,
                )
                # if the current block is less than 3 hours and the ratio to previous and post blocks is less than 80%
                if current < 3 and current / (prev + post) < 0.8:
                    df["wear"][
                        df.block == blocks[i][0]
                    ] = 0  # rescore the wear period as non wear
                # if the current block is less than 6 hours and the ratio to previous and post blocks is less than 30%
                elif current < 6 and current / (prev + post) < 0.3:
                    df["wear"][
                        df.block == blocks[i][0]
                    ] = 0  # rescore the wear period as non wear
        df.drop(columns=["block"], inplace=True)
        return df

    def rescore_last_day(self, df):
        """
        Rescores wear detection data saved in a pandas dataframe. Rescoring rules specific to the
        last day of recorded data.

        :param df: pandas dataframe containing wear/nonwear predictions
        :return: rescored pandas dataframe of wear/nonwear predictions
        """
        # group classifications into wear and nonwear blocks
        df["block"] = (df.wear.diff().ne(0)).cumsum()
        blocks = list(df.groupby("block"))

        # get the start index of the last day
        last_day_index = df.index[-1] - pd.to_timedelta("24h")

        # iterate through blocks
        for i in range(1, len(blocks)):
            wear = blocks[i][1]["wear"].values[
                0
            ]  # get whether or not the block is wear
            if (
                wear and blocks[i][1].index[0] > last_day_index
            ):  # if wear, and it's the last day
                # get hour lengths of the previous and current blocks
                prev, current = len(blocks[i - 1][1]) * 0.25, len(blocks[i][1]) * 0.25
                # if the current block is less than 3 hours and the previous block is greater or equal to 1 hour
                if current < 3 and prev >= 1:
                    df["wear"][
                        df.block == blocks[i][0]
                    ] = 0  # rescore the wear period as non wear
        df.drop(columns=["block"], inplace=True)
        return df


class ColeKripke:
    """
    Runs sleep wake detection on epoch level activity data. Epochs are 1 minute long and activity is represented
    by an activity index.
    """

    def __init__(self, activity_index):
        """
        Initialization of the class

        :param activity_index: pandas dataframe of epoch level activity index values
        """
        self.activity_index = activity_index
        self.predictions = None

    def predict(self, sf=np.array(0.193125)):
        """
        Runs the prediction of sleep wake states based on activity index data.

        :param sf: scale factor to use for the predictions (default corresponds to scale factor optimized for use with
        the activity index, if other activity measures are desired the scale factor can be modified or optimized.)
        The recommended range for the scale factor is between 0.1 and 0.25 depending on the sensitivity to activity
        desired, and possibly the population being observed.

        :return: rescored predictions
        """
        kernel = (
            sf
            * np.array([4.64, 6.87, 3.75, 5.07, 16.19, 5.84, 4.024, 0.00, 0.00])[::-1]
        )
        scores = np.convolve(self.activity_index, kernel, "same")
        scores[scores >= 0.5] = 1
        scores[scores < 0.5] = 0

        # rescore the original predictions
        self.rescore(scores)
        return self.predictions

    def rescore(self, predictions):
        """
        Application of Webster's rescoring rules as described in the Cole-Kripke paper.

        :param predictions: array of predictions
        :return: rescored predictions
        """
        rescored = predictions.copy()
        # rules a through c
        wake_bin = 0
        for t in range(len(rescored)):
            if rescored[t] == 1:
                wake_bin += 1
            else:
                if (
                    14 < wake_bin
                ):  # rule c: at least 15 minutes of wake, next 4 minutes of sleep get rescored
                    rescored[t : t + 4] = 1.0
                elif (
                    9 < wake_bin < 15
                ):  # rule b: at least 10 minutes of wake, next 3 minutes of sleep get rescored
                    rescored[t : t + 3] = 1.0
                elif (
                    3 < wake_bin < 10
                ):  # rule a: at least 4 minutes of wake, next 1 minute of sleep gets rescored
                    rescored[t] = 1.0
                wake_bin = 0
        # rule d: 6 minutes or less of sleep surrounded by at least 10 minutes of wake on each side gets rescored
        sleep_bin = 0
        start_ind = 0
        for t in range(10, len(rescored) - 10):
            if rescored[t] == 0:
                sleep_bin += 1
                if sleep_bin == 1:
                    start_ind = t
            else:
                if 0 < sleep_bin <= 6:
                    if (
                        sum(rescored[start_ind - 10 : start_ind]) == 10.0
                        and sum(rescored[t : t + 10]) == 10.0
                    ):
                        rescored[start_ind:t] = 1.0
                sleep_bin = 0
        self.predictions = rescored


def band_pass_filter(
    data_df, sampling_rate, bp_cutoff, order, channels=["X", "Y", "Z"]
):
    """
    Band-pass filter a given sensor signal.

    :param data_df: dataframe housing sensor signals
    :param sampling_rate: sampling rate of signal
    :param bp_cutoff: filter cutoffs
    :param order: filter order
    :param channels: channels of signal to filter
    :return: dataframe of raw and filtered data
    """
    data = data_df[channels].values

    # Calculate the critical frequency (radians/sample) based on cutoff frequency (Hz) and sampling rate (Hz)
    critical_frequency = [
        bp_cutoff[0] * 2.0 / sampling_rate,
        bp_cutoff[1] * 2.0 / sampling_rate,
    ]

    # Get the numerator (b) and denominator (a) of the IIR filter
    [b, a] = signal.butter(
        N=order, Wn=critical_frequency, btype="bandpass", analog=False
    )

    # Apply filter to raw data
    bp_filtered_data = signal.filtfilt(b, a, data, padlen=10, axis=0)

    new_channel_labels = [ax + "_bp_filt_" + str(bp_cutoff) for ax in channels]

    data_df[new_channel_labels] = pd.DataFrame(bp_filtered_data)

    return data_df


def activity_index(signal_df, channels=["X", "Y", "Z"]):
    """
    Compute activity index of sensor signals.

    :param signal_df: dataframe housing desired sensor signals
    :param channels: channels of signal to compute activity index
    :return: dataframe housing calculated activity index
    """
    ai_df = pd.DataFrame()
    ai_df["activity_index"] = [np.var(signal_df[channels], axis=0).mean() ** 0.5]
    return ai_df


def bin2df(full_path):
    """

    Reads geneactiv .bin files into a pandas dataframe

    :param full_path: full path to geneactiv .bin file

    :return decode: pandas dataframe of GA data

    """
    with open(full_path, "rb") as in_file:
        full_line = in_file.readline()
        count = 0
        fs = ""
        df = []
        while full_line:
            full_line = in_file.readline()
            line = full_line[:].split("\r\n")[0]
            count += 1
            if count < 60:
                if "x gain" in line:
                    x_gain = int(line.split(":")[-1])

                if "x offset" in line:
                    x_offset = int(line.split(":")[-1])

                if "y gain" in line:
                    y_gain = int(line.split(":")[-1])

                if "y offset" in line:
                    y_offset = int(line.split(":")[-1])

                if "z gain" in line:
                    z_gain = int(line.split(":")[-1])

                if "z offset" in line:
                    z_offset = int(line.split(":")[-1])

                if "Volts" in line:
                    volts = int(line.split(":")[-1])

                if "Lux" in line:
                    lux = int(line.split(":")[-1])

            if "Page Time:" in line:
                time = pd.to_datetime(
                    ":".join(line.split(":")[1:])[0:-2], format="%Y-%m-%d %H:%M:%S:%f"
                )

            if "Temperature:" in line:
                temp = float(line.split(":")[-1])

            if not fs:
                if "Measurement Frequency:" in line:
                    fs = float(line.split(":")[-1].split(" ")[0])
                    offset = np.array([1 / fs] * 300) * np.arange(0, 300)
                    delta = pd.to_timedelta(offset, unit="s")

            if len(line) == 3600:
                # hex to bin
                hexes = struct.unpack("12s " * 300, line)
                bins = (
                    struct.unpack(
                        "12s 12s 12s 10s 1s 1s", bin(int(hx, 16))[2:].zfill(48)
                    )
                    for hx in hexes
                )
                decode = pd.DataFrame(
                    bins,
                    columns=["X", "Y", "Z", "Light", "Button", "_"],
                    index=pd.DatetimeIndex([time] * 300) + delta,
                )

                # binary to decimal and calibration
                decode.X = decode.X.apply(
                    lambda x: round(
                        (BitArray(bin=x).int * 100.0 - x_offset) / x_gain, 4
                    )
                )
                decode.Y = decode.Y.apply(
                    lambda x: round(
                        (BitArray(bin=x).int * 100.0 - y_offset) / y_gain, 4
                    )
                )
                decode.Z = decode.Z.apply(
                    lambda x: round(
                        (BitArray(bin=x).int * 100.0 - z_offset) / z_gain, 4
                    )
                )
                decode.Light = decode.Light.apply(lambda x: int(x, 2) * lux / volts)
                decode["T"] = temp
                df.append(decode)

        df = pd.concat(df, axis=0)
        df.index.name = "Time"
        return df[["X", "Y", "Z", "T"]]
