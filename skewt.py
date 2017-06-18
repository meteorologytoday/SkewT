from matplotlib import rc;
import matplotlib.pyplot as pplt;
import matplotlib.font_manager as fm;
import matplotlib as mplt;
import sys
import numpy as np

from thermodynamics import *

#matplotlib.use('Agg')

rc('font', **{'family':'sans-serif', 'serif': 'Bitstream Vera Serif', 'sans-serif': 'MS Reference Sans Serif', 'size': 20.0, 'weight' : 100})
rc('axes', **{'labelsize': 20.0, 'labelweight': 100})
rc('mathtext', **{'fontset':'stixsans'});
default_linewidth = 2.0;
default_ticksize = 10.0;
mplt.rcParams['lines.linewidth'] =   default_linewidth;
mplt.rcParams['axes.linewidth'] =    default_linewidth;
mplt.rcParams['xtick.major.size'] =  default_ticksize;
mplt.rcParams['xtick.major.width'] = default_linewidth;
mplt.rcParams['ytick.major.size'] =  default_ticksize;
mplt.rcParams['ytick.major.width'] = default_linewidth;
mplt.rcParams['xtick.major.pad']='12';
mplt.rcParams['ytick.major.pad']='12';
print("%s: %d"%(fm.FontProperties().get_name(),fm.FontProperties().get_weight()));


weight = -30.0


p_range = np.array([100e2, p_ground])
T_range = np.array([-35, 50]) + zeroK
logp_range = np.log(p_range)

p_vec = np.linspace(p_range[0], p_range[1], num=100)
logp_vec = np.log(p_vec)
wlogp_vec = weight * logp_vec

Tlogp_range = T_range + weight * np.log(p_ground)
print(T_range)
print(logp_range)
print(Tlogp_range)
p_lines = np.array([1000, 850, 700, 500, 400, 300, 250, 200, 150, 100]) * 100.0
T_lines = np.arange(-30, 50, 10) + zeroK 
dry_lines = np.arange(-30, 180, 10) + zeroK
wet_lines = np.arange(10, 180, 10) + zeroK

# calculate figure size
graph_size = [4.0, 6.0]
space = {
    'wspace': 2.0,
    'hspace': 1.8,
    'tspace': 1.0,
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
print(Tlogp_range)
fig = pplt.figure(figsize=figsize)
ax = fig.add_axes(rect, autoscale_on=False)

ax.set_yticks(np.log(p_lines))
ax.set_yticklabels(["%d" % int(p/100) for p in p_lines])

ax.set_xticks(T_lines + weight * logp_range[1])
ax.set_xticklabels(["%d" % np.round(T-zeroK) for T in T_lines])
ax.set_xlim(Tlogp_range)
ax.set_ylim(logp_range)
ax.invert_yaxis()

# p line
for p in p_lines:
	ax.plot(Tlogp_range, [np.log(p)]*2, color='k', linewidth=1)


# dry line:
for theta in dry_lines: 
	T_vec = (p_vec / p_ref) ** kappa * theta
	Tlogp_vec = T_vec + wlogp_vec
#	for i,_ in enumerate(T_vec):
#		print("p: %.2f, T: %.2f, theta: %.2f, Tlogp: %.2f or %.2f" % (p_vec[i], T_vec[i], theta, T_vec[i] + weight * np.log(p_vec[i]), Tlogp_vec[i]))#wlogp_vec[i]))
	ax.plot(Tlogp_vec, logp_vec, color='k', linewidth=1)

# wet line:
for theta_e in wet_lines: 
	T_vec = np.array([solve_T_given_theta_e_and_p(theta_e, p) for p in p_vec])
	Tlogp_vec = T_vec + wlogp_vec
	ax.plot(Tlogp_vec, logp_vec, color='k', linewidth=1, dashes=(10,5))

#ax2 = ax.twiny()
#ax2.set_xlim(np.amin(wlogp_vec), np.amax(wlogp_vec))
#ax2.set_xlim([-50, 50])
#ax2.plot(wlogp_vec, logp_vec)
#ax2.set_ylim([0,1])
#ax2.set_yticks([])

pplt.show()
