# JetCorePersistence
Library of several R and Python scripts to calculate eddy feedback, persistence and NAO / SNAO. Used in
`Understanding Persistent Summer North Atlantic Jet Fluctuations: The Role of Dynamical Feedbacks in Model and Observational Discrepancies` by Albert Ossó and Florian Ennemoser, Wegener Center for Climate and Global Change, University of Graz.

## FeedbackCalc
* NAO_cmip6.R: Calculates the first EOF of the NATL SLP using the functions eof.mca.R and cov4gappy.R
* SNAO_structures.R: Project the Vorticity equation terms onto the NAOI. Used to produce Fig.2
* Forcing_timeseries.R:  Calculates time series M, and Y
* Eddy_feedback.R: Calculates the feedback parameter strength. Used also to produce several plots explain in the file.


## Base
* preprocessing.py: Preprocesses netcdf (.nc) files in order to fit calculation setup.
* calculation.py: Implements functions to find latitudinal maxima, jet core latitudes, and jet latitudinal index (JLI), Eddy Kinetic Energy (EKE) and vectorized linear regression. Additionally function to calculate persistence based on latitudinal differences.
* postprocessing.py: plotting functionality
* parse_config.py: parses .toml files from Configuration folder
* main.py: calls functions from parse_config or individual functions from preprocessing/calculation and postprocessing.
