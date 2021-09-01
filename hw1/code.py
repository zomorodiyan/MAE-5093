import numpy as np
import matplotlib.pyplot as plt

def lagrange(x,points):
	#points.sort(key=lambda x: x[0])
	sum = 0
	for j in range(len(points)):
		tmp =1
		for i in range(len(points)):
			if(i!=j):
				tmp *= (x-points[i,0])/(points[j,0]-points[i,0]) 	
		sum += points[j,1]*tmp
	return sum
points = np.array([[1.0,0.038],[0.95,0.042],[0.81,0.058],[0.59,0.104],[0.31,0.295],[-1.0,0.038],[-0.95,0.042],[-0.81,0.058],[-0.59,0.104],[-0.31,0.295],[0.0,1.0]])
x = 0.7
#points = np.array([[0.0,0.0],[1.0,1.0],[3.0,2.0]])

#print('question 1 test! it should be -0.226 and it is: ', lagrange(x,points))


years_tuition = np.array([[1998,21300],[1999,23057],[2000,24441],[2001,25917],[2002,27204],[2003,28564],[2004,29847],[2005,31200],[2006,32994],[2007,34800],[2008,36030]])
plt.style.use('seaborn-whitegrid')
plt.plot(years_tuition[0,:], years_tuition[1,:], 'o', color='black');
