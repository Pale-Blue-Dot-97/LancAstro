
import random
import Dependencies.DataLoad as dl
import Plot2D as laplt
import MultiFig as mf

DATALABELS = []
POINTSTYLES = []

band_nums = ["427", "464", "484", "505", "527", "574", "624", "679", "709", "711", "738", "767"]
path = ["Plot_data\\"]
filenames = []
columnnames = ["Log_Mass_Center", "Log_Phi_no_completeness", "Log_Phi_Error_no_completeness_minus",
               "Log_Phi_Error_no_completeness_plus"]

for i in band_nums:
    filenames.append("Plot_data_Band%s.fits" % i)

columns = []

for i in range(len(filenames)):
    columns.append(columnnames)

data = dl.data_load(filenames, columns, path)

x = []
y = []
y_err = []

for i in range(len(filenames)):
    x.append([data[4 * i]])
    y.append([data[4 * i + 1]])

for i in range(len(band_nums)):
    y_err.append([[data[4 * i + 2], data[4 * i + 3]]])

pointstyle = ['o']
colours = ['k']

for i in range(len(band_nums)):
    r = lambda: random.randint(0, 255)
    random_colour = '#%02X%02X%02X' % (r(), r(), r())
    if i < len(colours):
        colours[i] = [random_colour]
    else:
        colours.append([random_colour])

for i in range(len(band_nums)):
    POINTSTYLES.append(pointstyle)
    DATALABELS.append([band_nums[i]])

shape = [[1, 2, 3],
         [4, 5, 6],
         [7, 8, 9],
         [10, 11, 12]]

mf.create_grid(x, y, shape, y_error=y_err, DATALABELS=DATALABELS, POINTSTYLES=POINTSTYLES, COLOURS=colours,
               figsize=(18, 25), x_label="Log_Mass_Center", y_label="Log_Phi")
