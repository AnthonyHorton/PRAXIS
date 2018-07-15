
# coding: utf-8

# In[1]:

import pylab as plt
import numpy as np
from matplotlib.patches import RegularPolygon


# In[ ]:


#hexArray[idx] = [centreX,centreY,flatRad,flux]
#fibre number = idx+1
hexArray = np.array([[]])
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]
hexArray[0] = [0,0,0.055,0]


# In[35]:

get_ipython().magic(u'')


# In[20]:

fig, ax = plt.subplots(1)
ax.set_aspect('equal')

for x, y, c, l in zip([0,0], [0,1], ['0.0','1.0'], ['red','blue']):
#     color = c[0]  # matplotlib understands lower case words for colours
    hex = RegularPolygon((x, y), numVertices=6, radius=2. / 3., 
                         orientation=np.radians(30), 
                         facecolor=c,  edgecolor='k')
    ax.add_patch(hex)
    # Also add a text label
    ax.text(x, y+0.2, l[0], ha='center', va='center', size=20)
plt.show()


# In[18]:

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(0)
n = 100
x = np.random.standard_normal(n)
y = np.random.standard_normal(n)
xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()

fig, ax = plt.subplots()
# fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)
hb = ax.hexbin(x, y, gridsize=500, cmap='inferno')
ax.axis([xmin, xmax, ymin, ymax])
ax.set_title("Hexagon binning")
# cb = fig.colorbar(hb, ax=ax)
# cb.set_label('counts')


plt.show()


# In[43]:

x


# In[44]:

x.shape


# In[47]:

plt.plot(y)

plt.show()


# In[ ]:



