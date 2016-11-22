n_elements = 10
twist = ['t'+str(n) for n in xrange(n_elements)]
chord = ['c'+str(n) for n in xrange(n_elements)]
omega = ['omega']
x = omega + twist + chord
g = ['0']*(2*(n_elements-1)+1)
g[0] = 'omega_con'
g_indx = 1
for i in xrange(len(twist)-1):
    g[g_indx] = twist[i] + '>' + twist[i+1]
    g_indx += 1
for i in xrange(len(chord)-1):
    g[g_indx] = chord[i] + '>' + chord[i+1]
    g_indx += 1

print g