import sys, getopt
import numpy as np
from Sounding import Sounding
from thermodynamics import *
import os.path
import matplotlib as mplt

def usage():
	print("Usage:")

try:
	opts, args = getopt.getopt(sys.argv[1:], "f:o:t:h", ["output-file=", "input-file=", "title=", "help"])
except getopt.GetoptError as err:
	# print help information and exit:
	print(str(err))  # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

title = ""
input_fname = None
output_fname = None
for o, a in opts:
	if o in ("-f", "--input-file"):
		input_fname = a
	elif o in ("-o", "--output-file"):
		output_fname = a
	elif o in ("-h", "--help"):
		usage()
		sys.exit()
	elif o in ("-t", "--title"):
		title = a 
	else:
		assert False, "unhandled option"


if output_fname is not None:
	mplt.use('Agg')
	print("Will output file: %s" % (output_fname,))

from matplotlib import rc
import matplotlib.pyplot as pplt
import matplotlib.font_manager as fm
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
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
#print("%s: %d"%(fm.FontProperties().get_name(),fm.FontProperties().get_weight()));

sounding = Sounding(input_fname) if input_fname is not None else None

weight = -35.44
T_color = '#dddddd'
p_ground = p_ref
p_range = np.array([100e2, 1100e2])
T_range = np.array([-35, 55]) + zeroK
logp_range = np.log(p_range)

p_vec  = np.linspace(p_range[0], p_ground, num=1000)
p2_vec = np.linspace(p_range[0], p_range[1], num=1000)
logp_vec = np.log(p_vec)
wlogp_vec = weight * logp_vec

Tlogp_range = T_range + weight * np.log(p_ground)
cape_range = np.array([-500, 3000.0])
barb_x = Tlogp_range[-1]+8

cape_ticks = np.array([-500, 0, 1000, 2000, 3000])
p_lines = np.array([1000, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100]) * 100.0
T_lines = np.arange(-100, 60, 10) + zeroK 
dry_lines = np.arange(-30, 180, 10) + zeroK
wet_lines = np.arange(10, 180, 10) + zeroK
mixing_ratio_lines = np.array([0.7, 1, 2, 4, 7, 10, 16, 24, 32, 40]) / 1000.0

# calculate figure size
main_graph_size = np.array([9.0, 12.0])
cape_graph_size = np.array([3.0, 12.0])
space = {
    'wspace': 2.0,
    'hspace': 1.8,
    'tspace': 2.0,
    'bspace': 2.0,
    'lspace': 2.0,
    'rspace': 1.5,
    'break' : 0.1
}

figsize = [ main_graph_size[0] + cape_graph_size[0] + space['wspace'] + space['lspace'] + space['rspace'],
            main_graph_size[1] + space['tspace'] + space['bspace'] ]

# create main axes
main_rect = (
    space['lspace'] / figsize[0],
    space['bspace'] / figsize[1],
    main_graph_size[0]   / figsize[0],
    main_graph_size[1]   / figsize[1]
)

cape_rect = (
    main_rect[0] + main_rect[2] + space['hspace'] / figsize[0],
    main_rect[1], 
	cape_graph_size[0]   / figsize[0],
    cape_graph_size[1]   / figsize[1]
)
print("Creating canvas...")
fig = pplt.figure(figsize=figsize)
ax      = fig.add_axes(main_rect, autoscale_on=False)
skew_ax = fig.add_axes(cape_rect, autoscale_on=False)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

hshrink = [0.7, 0.5]
wshrink = 0.7
height = abs(logp_range[1] - logp_range[0])
width  = abs(Tlogp_range[1] - Tlogp_range[0])
pts = np.array([
	[Tlogp_range[1], logp_range[0]],
	[Tlogp_range[1], logp_range[0] + height * hshrink[0]],
	[Tlogp_range[0] + wshrink * width, logp_range[0] + hshrink[1] *height],
	[Tlogp_range[0] + wshrink * width, logp_range[0]],
])
ax.add_patch(Polygon(pts, facecolor='#ffffff', edgecolor='#ffffff', zorder=100)) 
ax.plot([Tlogp_range[0], pts[3][0], pts[2][0], pts[1][0], Tlogp_range[1]], [logp_range[0], logp_range[0], pts[2][1], pts[1][1], logp_range[1]], color='k', linewidth=default_linewidth, clip_on=False, zorder=101)


ax.set_yticks(np.log(p_lines))
ax.set_yticklabels(["%d" % int(p/100) for p in p_lines])

ax.set_xticks(T_lines + weight * logp_range[1])
ax.set_xticklabels(["%d" % np.round(T-zeroK) for T in T_lines])

ax.tick_params('both', length=0)

ax.set_xlim(Tlogp_range)
ax.set_ylim(logp_range)

ax.invert_yaxis()

skew_ax.set_xlim(cape_range)
skew_ax.tick_params('y', length=0)
skew_ax.tick_params('x', length=10)
skew_ax.set_xticks(cape_ticks)
skew_ax.set_xticklabels(['-0.5', '0', '1', '2', '3'])


skew_ax.plot([0,0], [logp_range[0], logp_range[1]], color='#aaaaaa', linewidth=1)
skew_ax.set_ylim(logp_range)
skew_ax.set_yticks(np.log(p_lines))
skew_ax.set_yticklabels(["%d" % int(p/100) for p in p_lines])

skew_ax.invert_yaxis()

# labels
ax.text(0.5, 1.02, title, verticalalignment='bottom', horizontalalignment='center', fontsize=25, transform=ax.transAxes)
ax.text(0.5, -0.08, r'Temperature [$\degree  \mathrm{C}$]', verticalalignment='top', horizontalalignment='center', fontsize=20, transform=ax.transAxes)
ax.text(-0.1, 0.5, r'Pressure [$\mathrm{hPa}$]', verticalalignment='center', horizontalalignment='right', fontsize=20, transform=ax.transAxes, rotation=90)

skew_ax.text(0.5, -0.08, r'CAPE integral [$\mathrm{kJ} \,\, \mathrm{kg}^{-1}$]', verticalalignment='top', horizontalalignment='center', fontsize=20, transform=skew_ax.transAxes)

# p line
for p in p_lines:
	ax.plot(Tlogp_range, [np.log(p)]*2, color='k', linewidth=1)
	skew_ax.plot(cape_range, [np.log(p)]*2, color='#aaaaaa', linewidth=1)

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
#for theta_e in wet_lines: 
#	T_vec = np.array([solve_T_given_theta_es_and_p(theta_e, p) for p in p_vec])
#	Tlogp_vec = T_vec + wlogp_vec
#	ax.plot(Tlogp_vec, logp_vec, color='k', linewidth=1, dashes=(10,5))

# wet line:
rev_p_vec = p_vec[::-1]
for theta_e in np.arange(0, 41, 5) + zeroK: 
	T_vec = gen_pseudo_adiabatic_line(theta_e, rev_p_vec)[::-1]
	Tlogp_vec = T_vec + wlogp_vec
	ax.plot(Tlogp_vec, logp_vec, color='k', linewidth=1, dashes=(10,5))

#
# mixing ratio line:
for mix_r in mixing_ratio_lines: 
	T_vec = np.array([inv_saturated_vapor_mass(mix_r, p) for p in p_vec])
	Tlogp_vec = T_vec + wlogp_vec
	ax.plot(Tlogp_vec, logp_vec, color='g', linewidth=1, dashes=(3,3))

	mix_r *= 1000
	ax.text(Tlogp_vec[-1], logp_vec[-1] + 0.03, ("%.1f" if mix_r < 1 else "%.0f") % (mix_r, ), color='g', horizontalalignment='center', verticalalignment='top', fontsize=10)

ax.text(zeroK + 42.0 + weight * np.log(p_ground), logp_vec[-1] + 0.03, r'g / kg', color='g', horizontalalignment='center', verticalalignment='top', fontsize=10)

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

	s_cape_int       = data['CAPE_int']
	s_cape_int_p10   = data['CAPE_int_p10']

	# parcel
	s_parcel_Tlogp = data['PARCEL_T'] + s_wlogp
	#s_parcel_Tlogp_p10 = data['PARCEL_T_p10'] + s_wlogp
	
	ax.plot(s_Tlogp, s_logp, color='b', linewidth=2)
	ax.plot(s_dewlogp, s_logp, color='r', linewidth=2)
	ax.plot(s_parcel_Tlogp, s_logp, color='#ee00ee', linewidth=2)
	#ax.plot(s_parcel_Tlogp_p10, s_logp, color='#ee00ee', linewidth=2)

	skew_ax.plot(s_cape_int    , s_logp, color='k')
	#skew_ax.plot(s_cape_int_p10, s_logp, color='r')

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

	beg = logp_range[1] - height
	if 'LCL' in data:
		ax.text(Tlogp_range[0] + width, beg, "LCL : %.0f" % (data['LCL']/100.0), horizontalalignment='right', verticalalignment='top', zorder=200)

	beg += height / 20.0
	if 'CAPE_int' in data:
		ax.text(Tlogp_range[0] + width, beg, "CAPE: %.0f" % (np.amax(data['CAPE_int']),), horizontalalignment='right', verticalalignment='top', zorder=200)

	beg += height / 20.0
	if 'CAPE_int_p10' in data:
		ax.text(Tlogp_range[0] + width, beg, "CAPEp10 Max: %.0f" % (np.amax(data['CAPE_int_p10']),), horizontalalignment='right', verticalalignment='top', zorder=200)
		beg += height / 20.0

if output_fname is not None:
	print("Output file: %s" % (output_fname,))
	fig.savefig(output_fname, dpi=300)

else:
	pplt.show()


