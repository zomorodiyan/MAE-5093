import numpy as np
import matplotlib.pyplot as plt

def lagrange(x,points):
    #points.sort(key=lambda x: x[0])
    sum = 0.0
    for j in range(len(points)):
        tmp = 1.0
        for i in range(len(points)):
            if i!=j:
                tmp *= (x-points[i,0])/(points[j,0]-points[i,0])
        sum += points[j,1]*tmp
    return sum
#points = np.array([[0.0,0.0],[1.0,1.0],[3.0,2.0]])
#print('question 1 test! it should be -0.226 and it is: ', lagrange(x,points))

def main():
    points1 = np.empty((11,2))
    points1[:,0] = np.array(\
            [-1.0,-0.95,-0.81,-0.59,-0.31,0,1.0,0.95,0.81,0.59,0.31])
    points1[:,1] = np.power(25*np.power(points1[:,0], 2)+1, -1)
    np.set_printoptions(precision=4)
    print('question 1 (a) lagrange(0.9, points) = ', format(lagrange(0.9,points1),".4f"))

    points2 = np.empty((21,2))
    points2[:,0] = np.arange(-1,1.0001,0.1)
    points2[:,1] = np.power(25*np.power(points2[:,0], 2)+1, -1)


    x1 = np.arange(-1,1.00001,0.01)
    y1 = np.empty(x1.shape[0])
    for i in range(x1.shape[0]):
        y1[i] = lagrange(x1[i],points1)
    x2 = np.arange(-1,1.00001,0.01)
    y2 = np.empty(x2.shape[0])
    for i in range(x2.shape[0]):
        y2[i] = lagrange(x2[i],points2)
    fig1 = plt.figure()
    plt.style.use('seaborn-whitegrid')
    plt.plot(x1,y1,'r')
    plt.plot(x2,y2,'b')

    axes = plt.gca()
    axes.set_ylim([-5,5])

    fig1.savefig('1_b')
    print('question 1 (b) a figure is saved in the code directory')

main()
