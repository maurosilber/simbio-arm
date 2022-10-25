import matplotlib.pyplot as plt
import numpy as np
from simbio.simulator import Simulator

from .. import ARM


t = np.linspace(0, 1, 10)
sim = Simulator(ARM)
sim.run(t).filter(like="monomer").plot()
plt.show()
