import pandas as pd
import struct
import numpy as np
from more_itertools import run_length
from bitstring import BitArray
from scipy import signal


def bin2df(full_path):
    """
    Reads geneactiv .bin files into a pandas dataframe.

    Parameters
    ----------
    full_path : str
        Full path to the geneactiv .bin file.

    Returns
    -------
    decode : dataframe
        Dataframe of decoded geneactiv data.

    """
    with open(full_path, "rb") as in_file:
        full_line = in_file.readline()
        count = 0
        fs = ""
        df = []
        while full_line:
            full_line = in_file.readline()
            line = full_line[:].split(b"\r\n")[0]
            count += 1
            if count < 60:
                if b"x gain" in line:
                    x_gain = int(line.split(b":")[-1])

                if b"x offset" in line:
                    x_offset = int(line.split(b":")[-1])

                if b"y gain" in line:
                    y_gain = int(line.split(b":")[-1])

                if b"y offset" in line:
                    y_offset = int(line.split(b":")[-1])

                if b"z gain" in line:
                    z_gain = int(line.split(b":")[-1])

                if b"z offset" in line:
                    z_offset = int(line.split(b":")[-1])

                if b"Volts" in line:
                    volts = int(line.split(b":")[-1])

                if b"Lux" in line:
                    lux = int(line.split(b":")[-1])

            if b"Page Time:" in line:
                time = pd.to_datetime(
                    ":".join(line.decode().split(":")[1:])[0:-2],
                    format="%Y-%m-%d %H:%M:%S:%f",
                )

            if b"Temperature:" in line:
                temp = float(line.split(b":")[-1])

            if not fs:
                if b"Measurement Frequency:" in line:
                    fs = float(line.split(b":")[-1].split(b" ")[0])
                    offset = np.array([1 / fs] * 300) * np.arange(0, 300)
                    delta = pd.to_timedelta(offset, unit="s")

            if len(line) == 3600:
                # hex to bin
                hexes = struct.unpack(b"12s " * 300, line)
                bins = (
                    struct.unpack(
                        b"12s 12s 12s 10s 1s 1s",
                        bin(int(hx, 16))[2:].zfill(48).encode(),
                    )
                    for hx in hexes
                )
                decode = pd.DataFrame(
                    bins,
                    columns=["X", "Y", "Z", "LUX", "Button", "_"],
                    index=pd.DatetimeIndex([time] * 300) + delta,
                )

                # binary to decimal and calibration
                decode.X = decode.X.apply(
                    lambda x: round(
                        (BitArray(bin=x.decode()).int * 100.0 - x_offset) / x_gain, 4
                    )
                )
                decode.Y = decode.Y.apply(
                    lambda x: round(
                        (BitArray(bin=x.decode()).int * 100.0 - y_offset) / y_gain, 4
                    )
                )
                decode.Z = decode.Z.apply(
                    lambda x: round(
                        (BitArray(bin=x.decode()).int * 100.0 - z_offset) / z_gain, 4
                    )
                )
                decode.LUX = decode.LUX.apply(lambda x: int(int(x, 2) * lux / volts))
                decode["T"] = temp
                df.append(decode)

        df = pd.concat(df, axis=0)
        df.index.name = "Time"
        return df[["X", "Y", "Z", "LUX", "T"]]


def non_overlapping_windows(df, window_size):
    """
    Slices input data frame into windows of specified size.

    Parameters
    ----------
    df : data frame
        Data frame to be windowed along the vertical axis.
    window_size : int
        Window size in number of samples.

    Returns
    -------
    windowed_data : numpy array
        Data windowed to the specified size.

    """
    data_mat = np.asarray(df)
    r, c = data_mat.shape
    num_windows = r // window_size
    data_mat = data_mat[0 : num_windows * window_size, :]
    windowed_data = data_mat.reshape(num_windows, window_size, c)
    return windowed_data


def activity_index(windowed_array):
    """
    Compute activity index of windowed tri-axis accelerometer signal.

    Parameters
    ----------
    windowed_array : array-like
        Full X,Y,Z tri-axis accelerometer signal with dimensions [number of windows, size of window, 3]

    Returns
    -------
    activity_indices : array-like
        Array of activity indices for input signal.

    """

    activity_indices = (np.var(windowed_array, axis=2).mean(axis=1) ** 0.5).reshape(
        -1, 1
    )
    return activity_indices


def high_pass_filter(windowed_array, sampling_rate, hp_cutoff, order):
    """
    High-pass filter a given windowed sensor signal.

    Parameters
    ----------
    windowed_array : array-like
        Sensor signal windowed with shape [number_of_windows, samples_per_window, number_of_data_channels].
    sampling_rate : float
    hp_cutoff : float
        High pass filter cutoff.
    order : int
        Order of the filter.

    Returns
    -------
    filtered : array-like
        High pass filtered data.

    """

    critical_frequency = [hp_cutoff * 2.0 / sampling_rate]

    # NUMERATOR AND DENOMINATOR OF IIR filter
    [b, a] = signal.butter(
        N=order, Wn=critical_frequency, btype="highpass", analog=False
    )

    # APPLY FILTER
    filtered = signal.filtfilt(b, a, windowed_array, padlen=10, axis=1)

    return filtered


def load_geneactiv_csv(full_path):
    """
    Loads geneactiv csv data into a pandas dataframe.

    Parameters
    ----------
    full_path : string
        Full path to a geneactiv accelerometer csv file.

    Returns
    -------
    df : data frame
        Pandas data frame of geneactiv accelerometer data including temperature and light.

    """
    df = pd.read_csv(
        full_path,
        index_col=0,
        skiprows=100,
        header=None,
        names=["Time", "X", "Y", "Z", "LUX", "Button", "T"],
        usecols=["Time", "X", "Y", "Z", "LUX", "T"],
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
    df.index = pd.to_datetime(df.index, format="%Y-%m-%d %H:%M:%S:%f").values
    return df[["X", "Y", "Z", "LUX", "T"]]


def downsample_by_mean(df, fs="50L"):
    """
    Downsamples a pandas data frame to a desired frequency by mean based aggregation.

    Parameters
    ----------
    df : data frame
        Data to be downsampled. Data frame must have pandas date time index.
    fs : str
        New sampling rate in string format.

    Returns
    -------
    downsampled data frame

    """
    return df.resample(fs).mean()


def downsample(df, factor):
    """
    Downsamples a pandas data frame to a desired frequency after using an antialiasing filter.

    Parameters
    ----------
    df : data frame
        Data to be downsampled.
    factor : int
        Factor by which to downsample.

    Returns
    -------
    df : data frame
        Downsampled data frame

    """

    # START AND STOP OF INDEX
    start = df.index[0]
    end = df.index[-1]

    # COLUMN NAMES
    cols = df.columns

    # DOWNSAMPLE
    ds = signal.decimate(df.values, q=factor, axis=0)

    # BUILD INDEX
    idx = pd.date_range(start, end, len(ds))

    # RECONSTRUCT DATA FRAME
    df = pd.DataFrame(ds)
    df.index = idx
    df.columns = cols
    return df


def make_dummy_data(freq="50L"):
    """
    Makes dummy sleep data. Default dummy data is 24 hours at 20hz. Sleep window from 10pm to 8am. Temperature set at 27
    degrees celsius.

     Parameters
    ----------
    freq : string
        Frequency expressed as string. 20hz Default expressed as "50L". 100hz as "10L".

    Returns
    -------
    input_df, start_sleep, end_sleep : 3-tuple of data frame, date-time object, date-time object
        Dataframe of dummy sleep data, start of sleep window, end of sleep window.

    """

    # START AND END DATES
    start = pd.to_datetime("2018-01-01 12:00:00:000", format="%Y-%m-%d %H:%M:%S:%f")
    end = pd.to_datetime("2018-01-02 12:00:00:000", format="%Y-%m-%d %H:%M:%S:%f")

    # SLEEP PERIOD
    start_sleep = pd.to_datetime(
        "2018-01-01 22:00:00:000", format="%Y-%m-%d %H:%M:%S:%f"
    )
    end_sleep = pd.to_datetime("2018-01-02 06:00:00:000", format="%Y-%m-%d %H:%M:%S:%f")

    # DUMMY DATE
    time = pd.date_range(start, end, freq=freq)
    data = np.random.uniform(-4, 5, [len(time), 5])

    # DATA FRAME
    input_df = pd.DataFrame(data)
    input_df.columns = ["X", "Y", "Z", "T", "LUX"]
    input_df.index = time

    # INSERT SLEEP PERIOD + TEMPERATURE VALUES
    input_df[start_sleep:end_sleep] = 0.001
    input_df["T"] = 27.0
    input_df["LUX"] = np.random.uniform(0, 80, [len(time)])

    return input_df, start_sleep, end_sleep


def total_sleep_time(predictions):
    """
    Calculates total sleep time on an array of sleep/wake predictions in one minute epochs.

    Parameters
    ----------
    predictions : array-like
        Binary sleep/wake predictions. Awake encoded as 1 and sleep as 0.

    Returns
    -------
    tst : int
        Total time spent asleep based on sleep/wake predictions provided.

    """
    tst = len(predictions) - predictions.sum()
    return int(tst)


def percent_time_asleep(predictions):
    """
    Calculates the percent of time spent asleep on an array of sleep/wake predictions in one minute epochs.

    Parameters
    ----------
    predictions : array-like
        Binary sleep/wake predictions. Awake encoded as 1 and sleep as 0.

    Returns
    -------
    pta : float
        Percentage of time spent asleep based on sleep/wake predictions provided.

    """
    pta = 100.0 * (len(predictions) - predictions.sum()) / float(len(predictions))
    return np.round(pta, decimals=3)


def number_of_wake_bouts(predictions):
    """
    Calculates the number of wake bouts present in an array of sleep/wake predictions in one minute epochs. Number of
    wake bouts is exclusive of first wake before sleep, and last wake after sleep.

    Parameters
    ----------
    predictions : array-like
        Binary sleep/wake predictions. Awake encoded as 1 and sleep as 0.

    Returns
    -------
    nwb : int
        Total number of wake bouts based on sleep/wake predictions provided.

    """
    first_sleep_epoch = predictions.argmin()
    last_sleep_epoch = predictions[::-1].argmin()
    predictions = predictions[first_sleep_epoch : -1 - last_sleep_epoch]
    bouts = list(run_length.encode(predictions))
    nwb = len([i for i in bouts if i[0] == 1])
    return nwb


def sleep_onset_latency(predictions):
    """
    Calculates sleep onset latency on an array of sleep/wake predictions in one minute epochs. This corresponds to the
    total number of minutes awake before the first sleep period.

    Parameters
    ----------
    predictions : array-like
        Binary sleep/wake predictions. Awake encoded as 1 and sleep as 0.

    Returns
    -------
    sol : int
        Total number of minutes spent awake before the first sleep period.

    """
    first_sleep_epoch = predictions.argmin()
    sol = predictions[0:first_sleep_epoch].sum()
    return int(sol)


def wake_after_sleep_onset(predictions):
    """
    Calculates number of minutes awake after first sleep period on an array of sleep/wake predictions in one minute
    epochs. Value is exclusive of the last wake period after sleep.

    Parameters
    ----------
    predictions : array-like
        Binary sleep/wake predictions. Awake encoded as 1 and sleep as 0.

    Returns
    -------
    waso : int
        Total number of minutes spent awake after the first sleep period.

    """
    first_sleep_epoch = predictions.argmin()
    last_sleep_epoch = predictions[::-1].argmin()
    waso = predictions[first_sleep_epoch : -1 - last_sleep_epoch]
    waso = waso.sum()

    return int(waso)
