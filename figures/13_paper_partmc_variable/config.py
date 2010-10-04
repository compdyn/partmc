import matplotlib
i_loop_max = 100

netcdf_dir = "../../scenarios/5_weighted/out"

figure_width_single = 8.4 / 2.54
figure_width_double = 16.45 / 2.54

i_weighting_schemes = 6

s_crit_1 = 0.01
s_crit_2 = 0.1
s_crit_3 = 0.3
s_crit_4 = 0.6

c_value = 0.95

matplotlib.rc('xtick', labelsize = 7)
matplotlib.rc('legend', fontsize = 7, borderpad = 0.7, borderaxespad = 1)
matplotlib.rc('font', size = 7, family = "serif",
              serif = ["Computer Modern Roman"])







partmc_tool_dir = "../../tool"

plume_data_dir = "/home/ching1/subversion/partmc/trunk/scenarios/3_condense/start"
parcel_data_dir = "/home/ching1/subversion/partmc/trunk/scenarios/3_condense/out"
parcel_regexp = "cond_08_%s_0001_([0-9]{8}).nc"
parcel_init_regexp = "cond_([0-9]{2})_%s_0001_00000001.nc"
parcel_series_regexp = "cond_%s_%s_0001_([0-9]{8}).nc"

output_dir = "../graphics"
data_dir = "data"

averaging_schemes_plume = [
    ("", "reference"),
    ("_comp", "composition"),
    ("_size", "diameter"),
    ("_both", "both"),
    ]

bc_num_range = {
    "": (1e3, 1e5),
    "_comp": (1e4, 1e6),
    "_size": (1e4, 1e6),
    "_both": (1e4, 1e7),
    }

scrit_num_range = {
    "": (1e2, 1e5),
    "_comp": (1e4, 1e6),
    "_size": (1e4, 1e6),
    "_both": (1e4, 1e7),
    }

averaging_schemes_parcel = [
    ("ref", "reference"),
    ("comp", "composition"),
    ("size", "diameter"),
    ("both", "both"),
    ]

plume_hours = [1, 7, 15, 48]
parcel_minutes = [
    {"mins": 0, "label_diam": 1e-1,
     "label_ha": "right", "label_va": "bottom"},
    {"mins": 2, "label_diam": 4e0,
     "label_ha": "right", "label_va": "bottom"},
    {"mins": 10, "label_diam": 13,
     "label_ha": "left", "label_va": "bottom"},
    ]

activated_diameter = 2e-6
