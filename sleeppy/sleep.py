import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.backends.backend_pdf
from sleeppy.majorrest import MajorRestPeriod
from sleeppy.colekripke import ColeKripke
from sleeppy.utils import (
    downsample,
    bin2df,
    high_pass_filter,
    activity_index,
    load_geneactiv_csv,
    non_overlapping_windows,
    total_sleep_time,
    percent_time_asleep,
    number_of_wake_bouts,
    wake_after_sleep_onset,
    sleep_onset_latency,
)

sns.set()


class SleepPy(object):
    """
    Processes raw GeneActiv accelerometer data from the wrist to determine various sleep metrics and endpoints.

    *** Support for data from other wrist worn accelerometers will be added in the future. To maximize generality an
    input format will be specified. As long as the input specifications are met, the package should produce the proper
    results. ***

    The sleep window detection of this package are based off of the following papers, as
    well as their implementation in the R package GGIR:

    van Hees V, Fang Z, Zhao J, Heywood J, Mirkes E, Sabia S, Migueles J (2019). GGIR: Raw Accelerometer Data Analysis.
    doi: 10.5281/zenodo.1051064, R package version 1.9-1, https://CRAN.R-project.org/package=GGIR.

    van Hees V, Fang Z, Langford J, Assah F, Mohammad Mirkes A, da Silva I, Trenell M, White T, Wareham N,
    Brage S (2014).
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

    """

    def __init__(
        self,
        input_file,
        results_directory,
        sampling_frequency,
        aws_object=None,
        start_buffer="0s",
        stop_buffer="0s",
        start_time="",
        stop_time="",
        temperature_threshold=25.0,
        minimum_rest_block=30,
        allowed_rest_break=60,
        minimum_rest_threshold=0.0,
        maximum_rest_threshold=1000.0,
        minimum_hours=6,
    ):
        """
        Class initialization.

        Parameters
        ----------
        input_file : string
            Full path to the file to be processed.
        results_directory : string
            Full path to the directory where the results should be saved.
        sampling_frequency : float
            Sampling frequency of the input data.
        aws_object : data object
            Data object to be processed from aws (in place of source file path).
        start_buffer : string
            Number of seconds to ignore from the beginning of the recording. format: "0s"
        stop_buffer : string
             Number of seconds to ignore at the end of the recording. format: "0s"
        start_time : string
            Time stamp from which to start looking at data. format: "%Y-%m-%d %H:%M:%S:%f"
        stop_time : string
            Time stamp at which to stop looking at data. format: "%Y-%m-%d %H:%M:%S:%f"
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
        minimum_hours : float
            Minimum number of hours required to consider a day useable.

        """
        if aws_object is not None:
            self.src = aws_object
        else:
            self.src = input_file
        self.extension = input_file.split(".")[-1]
        self.dst = results_directory
        self.src_name = input_file.split("/")[-1][0:-4]
        self.sub_dst = results_directory
        self.fs = sampling_frequency
        self.window_size = 60
        self.high_pass_cutoff = 0.25
        self.start_buffer = start_buffer
        self.stop_buffer = stop_buffer
        self.start_time = start_time
        self.stop_time = stop_time
        self.min_t = temperature_threshold
        self.minimum_rest_block = minimum_rest_block
        self.allowed_rest_break = allowed_rest_break
        self.minimum_rest_threshold = minimum_rest_threshold
        self.maximum_rest_threshold = maximum_rest_threshold
        self.minimum_hours = minimum_hours

        # DATA
        self.raw_days = []
        self.raw_days_to_plot = []
        self.arm_angle_days = []
        self.rest_period_days = []
        self.activity_index_days = []
        self.sleep_wake_prediction_days = []
        self.major_rest_periods = []
        self.endpoints = [
            [
                "total_sleep_time",
                "percent_time_asleep",
                "waso",
                "sleep_onset_latency",
                "number_wake_bouts",
            ]
        ]

    def run(self):
        """
        Runs the full package end to end on the input data.

        """
        print("Loading data...")
        self.load()
        print("Detecting major rest periods...")
        self.calculate_major_rest_periods()
        print("Making sleep predictions...")
        self.calculate_sleep_predictions()
        print("Calculating sleep endpoints...")
        self.calculate_endpoints()
        print("Building plots...")
        self.visualize()
        print("Exporting results...")
        self.export()

        return

    def load(self):
        """
        Loads full data into 24 hour days.

        """
        # FULL DATA
        df = (
            load_geneactiv_csv(self.src)
            if self.extension == "csv"
            else bin2df(self.src)
        )

        # DOWNSAMPLE
        if self.fs > 20.0:
            df = downsample(df, int(self.fs) // 20)
            self.fs = 20.0

        # APPLY BUFFERS
        df = df.loc[
            df.index[0]
            + pd.Timedelta(self.start_buffer) : df.index[-1]
            - pd.Timedelta(self.stop_buffer)
        ]

        # APPLY START AND STOP TIMES
        if self.start_time and self.stop_time:
            self.start_time = pd.to_datetime(
                self.start_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            self.stop_time = pd.to_datetime(
                self.stop_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            df = df.loc[self.start_time : self.stop_time]
        elif self.start_time:
            self.start_time = pd.to_datetime(
                self.start_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            df = df.loc[self.start_time :]
        elif self.stop_time:
            self.stop_time = pd.to_datetime(
                self.stop_time, format="%Y-%m-%d %H:%M:%S:%f"
            )
            df = df.loc[: self.stop_time]

        # SPLIT BY DAY
        self.raw_days = [
            day[1]
            for day in df.groupby(pd.Grouper(level=0, freq="24h", base=12))
            if ((len(day[1]) / float(self.fs)) / 3600.0) >= self.minimum_hours
        ]

        # DAYS FOR PLOTTING
        self.raw_days_to_plot = [day.resample("60s").median() for day in self.raw_days]
        return

    def calculate_major_rest_periods(self):
        """
        Gets the major rest periods for each 24 hour day.

        """
        for day in self.raw_days:
            mrp = MajorRestPeriod(day.copy(), self.fs)
            mrp.run()
            self.arm_angle_days.append(mrp.df_angle)
            self.rest_period_days.append(mrp.rest_periods)
            self.major_rest_periods.append(mrp.mrp)
        return

    def calculate_sleep_predictions(self):
        """
        Gets the sleep/wake predictions for each 24 hour day.

        """
        for day in self.raw_days:
            # PREPROCESS, ACTIVITY INDEX
            day = day.copy()[["X", "Y", "Z"]]
            windowed = non_overlapping_windows(day.copy(), int(60 * self.fs))
            filtered = high_pass_filter(
                windowed,
                sampling_rate=self.fs,
                hp_cutoff=self.high_pass_cutoff,
                order=3,
            )
            ai = activity_index(filtered)

            # SLEEP WAKE PREDICTIONS
            ck = ColeKripke(ai.flatten())
            predictions = ck.predict()

            # ADD TIME BACK
            idx = pd.date_range(day.index[0], day.index[-1], freq="60s")[
                0 : len(predictions)
            ]
            predictions = pd.DataFrame(predictions)
            predictions.index = idx
            ai = pd.DataFrame(ai)
            ai.index = idx

            # STORE DATA
            self.activity_index_days.append(ai)
            self.sleep_wake_prediction_days.append(predictions)
        return

    def calculate_endpoints(self):
        """
        Calculates the clinical sleep endpoints for each 24 hour day.

        """
        for day in range(len(self.sleep_wake_prediction_days)):
            mrp = self.major_rest_periods[day]
            if mrp:
                predictions = self.sleep_wake_prediction_days[day][mrp[0] : mrp[1]]
            else:
                predictions = self.sleep_wake_prediction_days[day]
            predictions = predictions.values
            self.endpoints.append(
                [
                    total_sleep_time(predictions),
                    percent_time_asleep(predictions),
                    wake_after_sleep_onset(predictions),
                    sleep_onset_latency(predictions),
                    number_of_wake_bouts(predictions),
                ]
            )
        return

    def visualize(self):
        """
        Builds all visual reports and exports them to pdf.

        """
        # OUTPUT PDF FOR ALL DAYS
        pdf = matplotlib.backends.backend_pdf.PdfPages(
            self.dst + "/{}_visual_report.pdf".format(self.src_name)
        )

        # PROCESS DAYS
        days = list(range(0, len(self.raw_days_to_plot)))
        for day in days:
            # LOAD DATA, GET SHARED AXIS
            raw = self.raw_days_to_plot[day]
            idx = pd.date_range(
                start=raw.index[0].replace(hour=12, minute=0, second=0, microsecond=0),
                periods=1440,
                freq="60s",
            )
            raw = raw.reindex(idx, fill_value=np.nan)

            # ACTIVITY INDEX
            act = self.activity_index_days[day]
            act = act.resample("60s").max()
            act = act.reindex(idx, fill_value=np.nan)

            # ARM ANGLE
            ang = self.arm_angle_days[day]
            ang = ang.resample("60s").max()
            ang = ang.reindex(idx, fill_value=np.nan)

            # SLEEP WAKE PREDICTIONS
            swp = self.sleep_wake_prediction_days[day]
            swp[swp == 0] = np.nan
            swp = swp.resample("60s").max()
            swp = swp.reindex(idx, fill_value=np.nan)

            # REST PERIODS
            rps = self.rest_period_days[day]
            rps[rps == 1] = np.nan
            rps[rps == 0] = 1
            rps = rps.resample("60s").max()
            rps = rps.reindex(idx, fill_value=np.nan)

            # DATAFRAME FOR STRAIGHT LINES
            df = swp.copy()
            df.columns = ["wake predictions"]
            df[["wake predictions"]] -= 0.25
            df["rest periods"] = rps.values - 0.75

            # ENDPOINTS FOR TABLE
            t_labels = (
                "Total Sleep Time(minutes)",
                "Percent Time Asleep",
                "Wake After Sleep Onset(minutes)",
                "Sleep Onset Latency(minutes)",
                "Number of Wake Bouts",
            )
            t_vals = [self.endpoints[day + 1]]

            # PLOTTING
            fig, (axt, ax0, ax1, ax2, ax3, ax4, ax5) = plt.subplots(
                7, 1, figsize=(30, 15)
            )
            plt.suptitle(
                "Visual Report for Source: {}\nDay: {}\nDate: {}".format(
                    self.src_name, day + 1, idx[0].date()
                ),
                fontsize=25,
            )
            hours = mdates.HourLocator(interval=1)
            h_fmt = mdates.DateFormatter("%H:%M")
            all_axes = (ax0, ax1, ax2, ax3, ax4, ax5)

            # PLOT TABLE
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

            # RAW DATA
            raw.rename(columns={"T": "Temperature", "LUX": "Light"}, inplace=True)
            raw[["X", "Y", "Z"]].plot(ax=ax0, lw=1).legend(
                bbox_to_anchor=(0, 1), fontsize=20
            )

            # TEMPERATURE
            raw[["Temperature"]].plot(
                ax=ax1, lw=1, color=sns.xkcd_rgb["pale red"]
            ).legend(bbox_to_anchor=(0, 1), fontsize=20)
            idx_first = raw.Temperature.first_valid_index()
            idx_last = raw.Temperature.last_valid_index()
            ax1.hlines(
                y=self.min_t,
                xmin=idx_first,
                xmax=idx_last,
                color="r",
                linestyle="--",
                lw=2,
            )
            props = dict(boxstyle="round", facecolor="lavender", alpha=0.35)
            textstr = u"max: {}\xb0C\nmin: {}\xb0C\nthresh: {}\xb0C".format(
                np.round(raw[["Temperature"]].max().values[0], decimals=2),
                np.round(raw[["Temperature"]].min().values[0], decimals=2),
                self.min_t,
            )
            ax1.text(
                0.005,
                0.95,
                textstr,
                transform=ax1.transAxes,
                fontsize=14,
                verticalalignment="top",
                bbox=props,
            )

            # LIGHT
            raw[["Light"]].plot(ax=ax2, lw=1, color=sns.xkcd_rgb["pale orange"]).legend(
                bbox_to_anchor=(0, 1), fontsize=20
            )

            # ACTIVITY INDEX
            act.plot(ax=ax3, lw=1, color="#6fc276").legend(
                labels=["activity"], bbox_to_anchor=(0, 0.75), fontsize=20
            )

            # ARM ANGLE
            ang.plot(ax=ax4, lw=1, color="#b36ff6").legend(
                labels=["arm angle"], bbox_to_anchor=(0, 0.75), fontsize=20
            )

            # SLEEP WAKE/REST PERIODS
            df.plot(ax=ax5, lw=10, x_compat=True, ylim=(0, 1)).legend(
                bbox_to_anchor=(0, 0.85), fontsize=20
            )

            # FORMAT
            plt.draw()
            count = 0
            for ax in all_axes:
                count += 1
                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.spines["bottom"].set_visible(False)
                ax.spines["left"].set_visible(False)
                ax.grid(False)
                if count < 6:
                    ax.get_xaxis().set_ticks([])
                ax.get_yaxis().set_ticks([])
            ax5.xaxis.set_major_locator(hours)
            ax5.xaxis.set_major_formatter(h_fmt)
            plt.subplots_adjust(wspace=0, hspace=0)
            fig.autofmt_xdate()
            for tick in ax5.xaxis.get_major_ticks():
                tick.label1.set_fontsize(16)
            pdf.savefig(fig.number)
        pdf.close()
        plt.close()

        # SUMMARY REPORT
        fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, 1, figsize=(12, 12))
        plt.suptitle("Summary Report for Source: {}".format(self.src_name), fontsize=16)
        all_axes = (ax0, ax1, ax2, ax3, ax4)

        # ENDPOINTS DATAFRAME
        endpoints = pd.DataFrame(self.endpoints[1:])
        endpoints.columns = self.endpoints[0]
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

        # TOTAL SLEEP TIME
        endpoints.total_sleep_time.plot.bar(
            ax=ax0, title="", color=["C" + str(i) for i in range(len(endpoints))]
        )
        # PERCENT TIME ASLEEP
        endpoints.percent_time_asleep.plot.bar(
            ax=ax1, title="", color=["C" + str(i) for i in range(len(endpoints))]
        )
        # WAKE AFTER SLEEP ONSET
        endpoints.waso.plot.bar(
            ax=ax2, title="", color=["C" + str(i) for i in range(len(endpoints))]
        )
        # SLEEP ONSET LATENCY
        endpoints.sleep_onset_latency.plot.bar(
            ax=ax3, title="", color=["C" + str(i) for i in range(len(endpoints))]
        )
        # NUMBER OF WAKE BOUTS
        endpoints.number_wake_bouts.plot.bar(
            ax=ax4, title="", color=["C" + str(i) for i in range(len(endpoints))]
        )
        # FORMAT
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
        plt.savefig(self.dst + "/{}_summary_report.pdf".format(self.src_name))
        plt.close()
        return

    def export(self):
        """
        Exports a results table to csv.

        """
        header = ["date"] + self.endpoints[0] + ["sleep_start", "sleep_end"]
        days = list(range(0, len(self.raw_days_to_plot)))
        export_df = []
        for day in days:
            export_df.append(
                [str(self.major_rest_periods[day][0].date())]
                + self.endpoints[day + 1]
                + [str(i) for i in self.major_rest_periods[day]]
            )
        export_df = pd.DataFrame(export_df)
        export_df.columns = header
        export_df.to_csv(
            self.dst + "/{}_output_table.csv".format(self.src_name), index=False
        )
        return
