import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from mpl_toolkits.mplot3d import Axes3D
import sys

fig = plt.figure() 
ax = fig.add_subplot(1,1,1, projection='3d')
ax.set_xlim3d([-5, 5])
ax.set_xlabel('X')
ax.set_ylim3d([-5, 5])
ax.set_ylabel('Y')
ax.set_zlim3d([-5, 5])
ax.set_zlabel('Z')
# ax.zscale("log")

plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

sc = ax.scatter([],[],[])
def init(): 
	return sc

# animation function 
def animate(i): 
	angle = i % 360
	xdata, ydata, zdata = [], [], []
	inp = sys.stdin.readline()
	if len(inp) == 0:
		anim.event_source.stop()

	inp = inp.rstrip()
	while inp != "":
		x,y,z = [ float(x) for x in inp.split() ][:-1][:3]
		xdata.append(x) 
		ydata.append(y) 
		zdata.append(z) 
		inp = sys.stdin.readline().rstrip()

	ax.view_init(30, int(angle))
	sc._offsets3d = (xdata, ydata, zdata)
	return sc 
	
# plt.axis('off') 
anim = animation.FuncAnimation(fig, animate, init_func=init, interval=50) 
plt.show()
