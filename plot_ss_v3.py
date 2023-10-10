##############################################################################
#                             PLOT SOLAR SYSTEM                              #
##############################################################################

'''
Written by Nicolas Kurtovic
'''

# Import
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astroquery.jplhorizons import Horizons
import datetime


##############################################################################
#                                  CLASS                                     #
##############################################################################


class SolarSystem:
    def __init__(self, initdate):
        '''
        initdate must be in the format 'YYYY-MM-DD'
        '''
        # List of planets
        self.planets  = []
        # Semi major axis
        self.semimaj = np.array([  0.3870321,  0.72326203, \
                                          1.,  1.52406417, \
                                   5.2038770,  9.57219251, \
                                  19.1644385, 30.18048128])
        # Periods in days
        self.periods  = np.array([  88,   225,   366,   687, \
                                  4331, 10747, 30589, 59800])
        # Number of days per step, and number of steps
        self.steps    = np.array([ 1,  1,  1,  1, \
                                   5, 10, 20, 30])
        self.nsteps   = (self.periods / self.steps).astype(int)
        # Set initial date
        self.initdate = initdate
        # Dictionary of positions
        self.pos      = {}

        # Iterate over the 8 planets of the Solay System
        for i in range(8):
            nasaid = i + 1
            # Create dictionary of positions for each planet
            self.pos[str(nasaid)] = {}
            # Find date for beginning of the orbit
            initorbit = (datetime.date.fromisoformat(initdate) - \
                         datetime.timedelta(days=int(self.periods[i]))).isoformat()
            # Get planet positions over its last orbit
            obj = Horizons(id=nasaid, \
                           epochs={'start': initorbit, \
                                   'stop':  self.initdate, \
                                   'step':  str(int(self.steps[i])) + 'd'}, \
                           id_type=None).vectors()
            # Obtain coordinates for plotting
            self.pos[str(nasaid)]['x'] = np.array([obj['x'][i] \
                                                  for i in range(self.nsteps[i])])
            self.pos[str(nasaid)]['y'] = np.array([obj['y'][i] \
                                                  for i in range(self.nsteps[i])])
            self.pos[str(nasaid)]['z'] = np.array([obj['z'][i] \
                                                  for i in range(self.nsteps[i])])
            # Close loop
            self.pos[str(nasaid)]['x'] = np.hstack([self.pos[str(nasaid)]['x'], \
                                                    self.pos[str(nasaid)]['x'][0]])
            # Close loop
            self.pos[str(nasaid)]['y'] = np.hstack([self.pos[str(nasaid)]['y'], \
                                                    self.pos[str(nasaid)]['y'][0]])
            # Close loop
            self.pos[str(nasaid)]['z'] = np.hstack([self.pos[str(nasaid)]['z'], \
                                                    self.pos[str(nasaid)]['z'][0]])


##############################################################################
#                                   MAIN                                     #
##############################################################################

#####################
# Set date
#####################

# Set an initial date in the format YYYY-MM-DD
initdate = '2023-10-10'   # October 10, 2023

# Calculate solar system
solars = SolarSystem(initdate)


#####################
# Ploting Parameters
#####################

# Geometry. Set to 0 for face on view. Choose any value between [0, 90]
inc = 55
# Power scaling of the planets semi-major axis. Use k=1 for linear scale
k = 0.2
# Figure horizontal size
figsize = 9


#####################
# Ploting Parameters
#####################

# Calculate scaling for the semi-major axis
scale = np.power(solars.semimaj, k)

# Create Figure
fig = plt.figure(figsize=[figsize, figsize * (2./3.)])
# Create ax
ax = plt.axes([0., 0., 1., 1.], \
              xlim=np.array([-31/scale[-1], 31/scale[-1]]), \
              ylim=np.array([-31/scale[-1], 31/scale[-1]]) * np.cos(np.deg2rad(inc)))
ax.set_aspect('equal')
ax.axis('off')

# Plot last orbit
for i in range(1, 9):
    j = str(i)
    ax.plot(solars.pos[j]['x'] / scale[i-1], \
            solars.pos[j]['y'] / scale[i-1] * np.cos(np.deg2rad(inc)) + \
            solars.pos[j]['z'] / scale[i-1] * np.sin(np.deg2rad(inc)), \
            '-k', lw=0.7, alpha=0.2)

# Plot last trace
for i in range(1, 9):
    j = str(i)
    last = np.min([int(len(solars.pos[j]['x']) / 10), 50])
    ax.plot(solars.pos[j]['x'][-last:-1] / scale[i-1], \
            solars.pos[j]['y'][-last:-1] / scale[i-1] * np.cos(np.deg2rad(inc)) + \
            solars.pos[j]['z'][-last:-1] / scale[i-1] * np.sin(np.deg2rad(inc)), \
            '--k', lw=0.9, alpha=0.7)

# Planet position on initdate
for i in range(1, 9):
    j = str(i)
    ax.plot(solars.pos[j]['x'][-1] / scale[i-1], \
            solars.pos[j]['y'][-1] / scale[i-1] * np.cos(np.deg2rad(inc)) + \
            solars.pos[j]['z'][-1] / scale[i-1] * np.sin(np.deg2rad(inc)), '.k')

fig.savefig('SolarSystem_'+initdate+'_clear.pdf')
plt.show()

