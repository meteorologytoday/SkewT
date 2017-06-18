import sys, getopt
import numpy as np


import matplotlib as mplt
#matplotlib.use('Agg')

from matplotlib import rc
import matplotlib.pyplot as pplt
import matplotlib.font_manager as fm

from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
from Sounding import Sounding
from thermodynamics import *

rc('font', **{'family':'sans-serif', 'serif': 'Bitstream Vera Serif', 'sans-serif': 'MS Reference Sans Serif', 'size': 15.0, 'weight' : 100})
rc('axes', **{'labelsize': 15.0, 'labelweight': 100})
rc('mathtext', **{'fontset':'stixsans'});
default_linewidth = 2.0;
default_ticksize = 10.0;
mplt.rcParams['lines.linewidth'] =   default_linewidth;
mplt.rcParams['axes.linewidth'] =	default_linewidth;
mplt.rcParams['xtick.major.size'] =  default_ticksize;
mplt.rcParams['xtick.major.width'] = default_linewidth;
mplt.rcParams['ytick.major.size'] =  default_ticksize;
mplt.rcParams['ytick.major.width'] = default_linewidth;
mplt.rcParams['xtick.major.pad']='12';
mplt.rcParams['ytick.major.pad']='12';
print("%s: %d"%(fm.FontProperties().get_name(),fm.FontProperties().get_weight()));

def usage():
	print("Usage:")

try:
	opts, args = getopt.getopt(sys.argv[1:], "f:t:h", ["file=", "title=", "help"])
except getopt.GetoptError as err:
	# print help information and exit:
	print(str(err))  # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

title = ""
fname = None
for o, a in opts:
	if o in ("-f", "--file"):
		fname = a
	elif o in ("-h", "--help"):
		usage()
		sys.exit()
	elif o in ("-t", "--title"):
		title = a 
	else:
		assert False, "unhandled option"

sounding = Sounding(fname) if fname is not None else None

weight = -30.0
T_color = '#dddddd'

p_range = np.array([100e2, p_ground])
T_range = np.array([-35, 50]) + zeroK
logp_range = np.log(p_range)

p_vec = np.linspace(p_range[0], p_range[1], num=100)
logp_vec = np.log(p_vec)
wlogp_vec = weight * logp_vec

Tlogp_range = T_range + weight * np.log(p_ground)
barb_x = Tlogp_range[-1]+8

p_lines = np.array([1000, 850, 700, 500, 400, 300, 250, 200, 150, 100]) * 100.0
T_lines = np.arange(-100, 60, 10) + zeroK 
dry_lines = np.arange(-30, 180, 10) + zeroK
wet_lines = np.arange(10, 180, 10) + zeroK
mixing_ratio_lines = np.array([0.7, 1, 2, 4, 7, 10, 16, 24, 32, 40]) / 1000.0
# calculate figure size
graph_size = np.array([4.0, 6.0]) * 2
space = {
    'wspace': 2.0,
    'hspace': 1.8,
    'tspace': 2.0,
    'bspace': 2.0,
    'lspace': 2.0,
    'rspace': 1.5,
    'break' : 0.1
}

figsize = [graph_size[0] + space['lspace'] + space['rspace'],
           graph_size[1] + space['tspace'] + space['bspace']]

# create main axes
rect = (
    space['lspace'] / figsize[0],
    space['bspace'] / figsize[1],
    graph_size[0]   / figsize[0],
    graph_size[1]   / figsize[1]
)
print("Creating canvas...")
fig = pplt.figure(figsize=figsize)
ax = fig.add_axes(rect, autoscale_on=False)

ax.set_yticks(np.log(p_lines))
ax.set_yticklabels(["%d" % int(p/100) for p in p_lines])

ax.set_xticks(T_lines + weight * logp_range[1])
ax.set_xticklabels(["%d" % np.round(T-zeroK) for T in T_lines])

ax.tick_params('both', length=0)

ax.set_xlim(Tlogp_range)
ax.set_ylim(logp_range)
ax.invert_yaxis()

# labels
ax.text(0.5, 1.02, title, verticalalignment='bottom', horizontalalignment='center', fontsize=25, transform=ax.transAxes)
ax.text(0.5, -0.08, r'$T - \mathrm{log} \, p$', verticalalignment='top', horizontalalignment='center', fontsize=20, transform=ax.transAxes)
ax.text(-0.1, 0.5, r'$\mathrm{log} \, p$', verticalalignment='center', horizontalalignment='right', fontsize=20, transform=ax.transAxes, rotation=90)


# p line
for p in p_lines:
	ax.plot(Tlogp_range, [np.log(p)]*2, color='k', linewidth=1)

# T line:
if len(T_lines) % 2 != 0:
	raise Error("There must be even numbers of T_lines.")

for i in range(int(len(T_lines)/2)):
	ii = 2*i
	T1, T2 = T_lines[ii], T_lines[ii+1]
	pts = np.array([
		[T1 + weight * logp_range[0], logp_range[0]],
		[T2 + weight * logp_range[0], logp_range[0]],
		[T2 + weight * logp_range[1], logp_range[1]],
		[T1 + weight * logp_range[1], logp_range[1]]
	])
	ax.add_patch(Polygon(pts, facecolor=T_color, edgecolor=T_color)) 

# dry line:
for theta in dry_lines: 
	T_vec = (p_vec / p_ref) ** kappa * theta
	Tlogp_vec = T_vec + wlogp_vec
	ax.plot(Tlogp_vec, logp_vec, color='k', linewidth=1)

# wet line:
for theta_e in wet_lines: 
	T_vec = np.array([solve_T_given_theta_e_and_p(theta_e, p) for p in p_vec])
	Tlogp_vec = T_vec + wlogp_vec
	ax.plot(Tlogp_vec, logp_vec, color='k', linewidth=1, dashes=(10,5))

# mixing ratio line:
for mix_r in mixing_ratio_lines: 
	T_vec = np.array([inv_saturated_vapor_mass(mix_r, p) for p in p_vec])
	Tlogp_vec = T_vec + wlogp_vec
	ax.plot(Tlogp_vec, logp_vec, color='g', linewidth=1, dashes=(3,3))

# wind barb line	
ax.plot([barb_x]*2, logp_range, color='k', linewidth=1, clip_on=False)

# Sounding part
if sounding is not None:
	data = sounding.data
	s_p = data['PRE']
	s_logp   = np.log(s_p)
	s_wlogp  = weight * s_logp
	
	s_Tlogp = data['TEMP']  + s_wlogp
	s_dewlogp = data['DEW'] + s_wlogp

	ax.plot(s_Tlogp, s_logp, color='b', linewidth=1)
	ax.plot(s_dewlogp, s_logp, color='r', linewidth=1)

	# wind
	upper_bound_i = 0
	for i, logp in enumerate(s_logp):
		if logp < logp_range[0]:
			upper_bound_i = i
			break	
	selected = slice(None,upper_bound_i,100)
	WS = data['WS'][selected]
	WD = data['WD'][selected]
	U = - WS * np.sin(WD * np.pi / 180.0)
	V = - WS * np.cos(WD * np.pi / 180.0)

	ax.barbs([barb_x]*len(U), s_logp[selected], U, V, WS, flagcolor='k', barbcolor='k', fill_empty=False, rounding=True, sizes=dict(emptybarb=0.1, spacing=0.1, height=0.3), clip_on=False)

fig.savefig("test.png", dpi=300)
pplt.show()

