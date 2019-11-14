import numpy as np


class ColeKripke(object):
    """
    Runs sleep wake detection on epoch level activity data. Epochs are 1 minute long and activity is represented
    by an activity index.

    """

    def __init__(self, activity_index, sf=0.243):

        """
        Class initialization.

        Parameters
        ----------
        activity_index : array-like
            Array of epoch level activity index values.
        sf : float
            Scale factor to use for the predictions (default corresponds to scale factor optimized for use with
            the activity index, if other activity measures are desired the scale factor can be modified or optimized.)
            The recommended range for the scale factor is between 0.15 and 0.3 depending on the sensitivity to activity
            desired, and possibly the population being observed.

        """

        self.activity_index = activity_index
        self.predictions = None
        self.rescored_predictions = None
        self.sf = sf

    def predict(self):
        """
        Runs full algorithm:
            Predict sleep/wake states using Cole-Kripke's optimal kernel with specified scale factor.
            Rescore predictions based on Webster's rules.

        Returns
        -------
        rescored_predictions : array-like

        """
        self.apply_kernel()
        return self.rescore()

    def apply_kernel(self):
        """
        Runs the prediction of sleep wake states based on activity index data.

        Returns
        -------
        predictions : array-like
            Sleep/wake predictions from a Cole-Kripke based algorithm.

        """

        kernel = (
            self.sf
            * np.array([4.64, 6.87, 3.75, 5.07, 16.19, 5.84, 4.024, 0.00, 0.00])[::-1]
        )
        scores = np.convolve(self.activity_index, kernel, "same")
        scores[scores >= 0.5] = 1
        scores[scores < 0.5] = 0
        self.predictions = scores
        return self.predictions

    def rescore(self):
        """
        Application of Webster's rescoring rules as described in the Cole-Kripke paper.

        Returns
        -------
        rescored : array-like
            Array of rescored sleep/wake predictions


        """

        rescored = self.predictions.copy()
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
        self.rescored_predictions = rescored
        return self.rescored_predictions
