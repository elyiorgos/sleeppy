---
title: 'SleepPy: A python package for sleep analysis from accelerometer data'
tags:
  - Python
  - SleepPy
  - Accelerometer
  - Wrist
  - Watch
  - GeneActiv
  - Digital
  - Wearables
authors:
  - name: Yiorgos Christakis
    orcid: 0000-0003-3371-8935
    affliation: 1
affiliations:
  - name: Pfizer, Inc.
  - index: 1
date: 12 July 2019
bibliography: paper.bib
---


# Introduction

Sleep quantity and quality are important measures in the evaluation of disease progression and regression across countless disease populations, as well as valuable measures of well-being in healthy populations. The use of wearable technologies to measure aspects of sleep quality and quantity in both healthy and clinical populations is still an evolving endeavor, and using wearable technologies, like wrist worn accelerometers, is increasingly becoming a valuable resource for extracting clinical endpoints related to sleep. 

``SleepPy`` is an open source python package that pulls together several published sleep analysis algorithms in a modular system to produce a suite of clinical endpoints to evaluate sleep. The package does the grunt work of processing raw accelerometer data from wrist-worn devices, running an array of analyses on the data and producing sleep reports on a 24-hour basis. The reports can be used to explore objective measures of sleep, gain insight into subject compliance(by highlighting obvious periods of recording time spent off-body), provide compelling visuals for the presentation of data to diverse audiences, and act as a visual debugging tool to evaluate algorithm performance, especially during development.

# Key Functionality
``SleepPy`` is designed for ease of use in a research environment, and therefore attempts to remove as much of the burden as possible from the user. It currently supports raw data collected with GeneActiv wrist-worn accelerometers, though support will be extended to a more general format for broader use. The package can be run with both GeneActiv .bin files, as well as the raw .csv outputs of the GeneActiv software. Processing the .bin files adds a non-trivial amount of time to each run of ``SleepPy``, but insofar as batch processing of .bin file extraction is not a feature of the GeneActiv software, the extra processing time is justified when dealing with a large number of files.

There are seven major steps that ``SleepPy`` goes through when processing data:

1. Splitting by day: Loads the input file and splits the resulting data into 24-hour segments, defined from noon to noon.

2. Activity index -`[@bai_di_xiao_evenson_lacroix_crainiceanu_buchner:2017]`: Calculates the activity index feature for 1 minute epochs for each of the days.

3. On-body detection -`[@hees_fang_zhao_heywood_mirkes_sabia_migueles:2019]`: Runs on-body/off-body detection for each of the days for debugging purposes

4. Major rest period -`[@hees_sabia_jones_wood_anderson_kivimäki_frayling_pack_bucan_trenell_et al.:2018]`: Calculates the major rest period (sleep window) for each day.

5. Sleep/wake -`[@cole_kripke_gruen_mullaney_gillin:1992]`: Runs sleep/wake predictions on each minute in the major rest period for each day.

6. Calculate endpoints: Using the above data, clinical endpoints for sleep are calculated for the major rest period of each day.

7. Visualize results: Creates a set of visual reports summarizing all findings for each day.

At the point of completion the user has access to all intermediate data, .csv files containing all relevant endpoint data, and two kinds of visual reports. The first report is a summary report of the endpoints across all days present in the input data. However, the most informative reports are the reports provided for each individual day, which resemble the image below (demo report):

![Demo Report](https://raw.githubusercontent.com/elyiorgos/sleeppy/master/sleeppy/demo/report_images/Visual_Results_Day_1.png?token=AECHRMACMGQTUEPUAPIXAEC5E53LG)

As shown, the report includes the source file name, the number of the day in the series of days provided, the start date of the data being shown, and a table of all calculated endpoints. Below the table is a graph of the data available during the 24-hour window specified. The subplots are set up to show the multiple forms that the data can take during the analysis. They are layed out as follows:

1. XYZ: The first level is simply the raw tri-axial accelerometer signal in X, Y, and Z.

2. Activity index: The second level is a plot of the minute by minute activity index values, which reflect the intensity of activity for each minute.

3. Arm-angle: The third level is a plot of the change in arm-angle over 24 hours, this is the feature used to determine the major rest period.

4. Wake: The fourth level represents sleep/wake predictions run on the entirety of the data, though only those predictions that fall under the major rest period are considered when computing clinical endpoints.

5. Rest periods: The fifth level is all periods of major rest detected by the major rest period algorithm. Only the largest, or major, rest period is used for calculating all clinical endpoints.

6. On-body: The sixth level is the raw predictions for on-body detection, without any filtering or rescoring.

7. On-body (rescored): The seventh level is the rescored predictions for on-body detection.

# Development
One of the driving factors in the design of ``SleepPy`` was the use of modularity to allow for the innovation and update of any algorithm or set of algorithms used in the package. The extraction of the activity index, for instance, can easily be modified to include the extraction of other features. Likewise, the sleep/wake detection algorithm can be replaced with another sleep/wake prediction module, as long as all input/output format criteria are maintained. The hope is that as new and better options become available, the package can be updated to provide the most reliable clinical endpoints possible. 

# Availability

The software is available as a pip installable package, as well as on GitHub at: <https://github.com/elyiorgos/sleeppy>

# Acknowledgements

This project was completed as part of research done for the Digital Medicine and Translational Imaging group at Pfizer, Inc. Specific guidance and assistance were provided throughout by Shyamal Patel and Nikhil Mahadevan, with input from the other members of the data science team including: Adan Rivas, Matthew Czech, Kara Chappie, Behnaz Rehzai.

# References

# Author
Yiorgos Christakis

# License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details