import matplotlib.pyplot as plt
import numpy as np

samples = (np.random.pareto((2.3-1), 1000) + 1)

count, bins, _ = plt.hist(samples, 100, density=True)
fit = (2.3-1)/ bins**(2.3)
plt.plot(bins, max(count)*fit/max(fit), linewidth=2, color='r')
plt.show()
