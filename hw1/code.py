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

def tdma(a,b,c,d):

    ac, bc, cc, dc = map(np.array, (a, b, c, d))
    ne = len(d)
    for i in range(1,ne):
        bc[i] = bc[i] - ac[i]*(cc[i-1]/bc[i-1])
        dc[i] = dc[i] - ac[i]*(dc[i-1]/bc[i-1])
    uc = bc
    uc[-1] = dc[-1]/bc[-1]

    for i in range(ne-2,-1,-1):
        uc[i] = (dc[i] - cc[i]*uc[i+1])/bc[i]

    del ac, bc, cc, dc

    return uc

def main():
#------------------------------question 1----------------------------------
    points1 = np.empty((11,2))
    points1[:,0] = np.array(\
            [-1.0,-0.95,-0.81,-0.59,-0.31,0,1.0,0.95,0.81,0.59,0.31])
    points1[:,1] = np.power(25*np.power(points1[:,0], 2)+1, -1)
    print('question 1 (a) lagrange(0.9, points) =', format(lagrange(0.9,points1),".4f"))

    points2 = np.empty((21,2))
    points2[:,0] = np.arange(-1,1.0001,0.1)
    points2[:,1] = np.power(25*np.power(points2[:,0], 2)+1, -1)

    x = np.arange(-1,1.00001,0.01)
    y1 = np.empty(x.shape[0]); y2 = np.empty(x.shape[0])
    for i in range(x.shape[0]):
        y1[i] = lagrange(x[i],points1)
        y2[i] = lagrange(x[i],points2)
    fig1 = plt.figure()
    plt.style.use('seaborn-whitegrid')
    plt.plot(x,y1,'r'); plt.plot(x,y2,'b')
    axes = plt.gca(); axes.set_ylim([-5,5])
    fig1.savefig('1_b')
    print('question 1 (b) figure saved in code directory')

#------------------------------question 7----------------------------------
    dat = np.array([[1998,21300],[1999,23057],[2000,24441]\
            ,[2001,25917],[2002,27204],[2003,28564],[2004,29847]\
            ,[2005,31200],[2006,32994],[2007,34800],[2008,36030]])
    fig2 = plt.figure()
    plt.style.use('seaborn-whitegrid')
    plt.plot(dat[:,0],dat[:,1])
    fig2.savefig('7_a')

#---------------------------j---question 8----------------------------------
    dat2 = np.array([[1993,12.0],[1995,12.7],[1997,13]\
            ,[1999,15.2],[2001,18.2],[2003,19.8],[2005,24.1]\
            ,[2007,28.1]])
    x = 2009

    #dat.sort(key=lambda x: x[0])

    delta=np.copy(dat[1:,0]-dat[0:-1,0])
    A = np.empty(dat.shape[0]); B = np.empty(dat.shape[0])
    C = np.empty(dat.shape[0]); D = np.empty(dat.shape[0])
    A[1:-1] = np.copy(delta[:-1]/6)
    B[1:-1] = np.copy((delta[:-1]+delta[1:])/3)
    C[1:-1] = np.copy(delta[1:]/6)
    D[1:-1] = np.copy((dat[2:,1]-dat[1:-1,1])/delta[1:]\
            -(dat[1:-1,1]-dat[:-2,1])/delta[:-1])
    A[0] = 0; A[-1] = 0
    B[0] = 1; B[-1] = 1
    C[0] = 0; C[-1] = 0
    D[0] = 0; D[-1] = 0

    G2nd = tdma(A,B,C,D)

    x = 2010
    interval = -1
    for i in range(dat.shape[0]):
        if x < dat[i,0]:
            interval = i
            break
    if interval == -1:
        interval = dat.shape[0]-2
    A = (dat[interval+1,0]-x)/(dat[interval+1,0]-dat[interval,0])
    B = 1 - A
    C = 1/6*(np.power(A,3)-A)*np.power((dat[interval+1,0]-dat[interval,0]),2)
    D = 1/6*(np.power(B,3)-B)*np.power((dat[interval+1,0]-dat[interval,0]),2)
    y = A*dat[interval,1] + B*dat[interval+1,1] + C*G2nd[interval]
    + D*G2nd[interval+1]
    print('y ',y)


main()
