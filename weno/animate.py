import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def animate(x, u):
    fig, ax = plt.subplots()
    line, = ax.plot(x, u[0], marker = 'o', linestyle = '-')

    def animate(i):
        line.set_ydata(u[i])  # update the data.
        return line,

    ani = animation.FuncAnimation(
        fig, animate, interval=100, blit=True, save_count=50)

    # To save the animation, use e.g.
    #
    #ani.save("RK3.mp4")
    #
    # or
    #
    #writer = animation.FFMpegWriter(
    #  fps=15, metadata=dict(artist='Me'), bitrate=1800)
    #ani.save("RK3.mp4", writer=writer)

    plt.show()

def main():
    x = np.linspace(0,10,101)
    u = np.zeros((100,101))
    for i in range(100):
        u[i]=np.sin(x+0.01*i)
    animate(x, u)

if __name__ == "__main__":
    main()
