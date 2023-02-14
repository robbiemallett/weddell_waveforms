# This repo analyses and plots KuKa radar waveforms

The code works from NetCDF files (nc files), which must be generated from the raw data using KuKaPy. The nc files should be placed in the data/kuka/required directory.

Because nc files are generated whenever the radar is running, there's a lot of data and we're only interested in a tiny bit of it. Specifically, we're only interested in data from each antenna when that antenna is pointing at the snow pit. The information about when that's happening was written in my notebook, and I've transcribed it into dictionaries within notebooks/pit_info.ipynb. Run this notebook first - it generates a pickle file of this information that you'll need later.

To make what I anticipate will be Figure 2, run notebooks/make_pit_plots.ipynb
