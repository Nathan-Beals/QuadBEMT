import numpy as np
import matplotlib.pyplot as plt

# data to plot
n_groups = 5
RPM_case1 = (5286, 5615, 6575, 6961, 5943)
power_case1 = (39.57, 39.51, 40.37, 42.14, 59.52)

RPM_case2 = (2496, 2601, 3184, 3543, 4740)
power_case2 = (38.12, 38.21, 38.42, 38.73, 65.85)

l_case1 = [r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.6}$',
           r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.4}$',
           r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.3}$',
           r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.2}$',
           r'$\mathrm{DA4002}$']

l_case2 = [r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.6}$',
           r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.4}$',
           r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.3}$',
           r'$\mathrm{(}\mathrm{c/R}\mathrm{)}_\mathrm{u}\mathrm{=}\mathrm{0.2}$',
           r'$\mathrm{Carroll}$']

# create plot
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
bar_width = 0.35
opacity = 0.8
index = np.arange(n_groups)+0.4*bar_width

rects1 = ax1.bar(index, RPM_case2, bar_width,
                 alpha=opacity,
                 color='b',
                 label='RPM')

rects2 = ax2.bar(index + bar_width, power_case2, bar_width,
                 alpha=opacity,
                 color='g',
                 label='power [W]')

ax1.set_xlabel("Case", fontsize=18)
ax1.set_ylabel('RPM', fontsize=18, color='b')
ax1.tick_params('y', colors='b')
ax2.set_ylabel("Power [W]", fontsize=18, color='g')
ax2.tick_params('y', colors='g')
plt.xticks(index + bar_width, l_case2, fontsize=18)
ax1.set_ylim([0, 8000])
ax2.set_ylim([0, 80])

plt.tight_layout()
plt.show()


# import numpy as np
# import matplotlib.pyplot as plt
#
# fig, ax1 = plt.subplots()
# t = np.arange(0.01, 10.0, 0.01)
# s1 = np.exp(t)
# ax1.plot(t, s1, 'b-')
# ax1.set_xlabel('time (s)')
# # Make the y-axis label, ticks and tick labels match the line color.
# ax1.set_ylabel('exp', color='b')
# ax1.tick_params('y', colors='b')
#
# ax2 = ax1.twinx()
# s2 = np.sin(2 * np.pi * t)
# ax2.plot(t, s2, 'r.')
# ax2.set_ylabel('sin', color='r')
# ax2.tick_params('y', colors='r')
#
# fig.tight_layout()
# plt.show()