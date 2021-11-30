import numpy as np
from animate import animate

# %%
def ux_coeff(start_term,n_terms):
    #start_term is the negative of the first term for u_i+0.5
    #start_term is the negative of the first term + 1 for u_i-0.5
    A = np.zeros((n_terms+1,n_terms+1))
    B = np.zeros((n_terms+1,1))
    B[1] = 1
    A[0] = np.ones((n_terms+1))

    for i in range(1,n_terms+1):
        for j in range(n_terms+1):
            A[i,j] = (-start_term-1+j)**i
    return np.linalg.inv(A)@B

# %%
def flux_coeff(start_term,n_terms):
    A = np.zeros((n_terms,n_terms))
    B = ux_coeff(start_term,n_terms)
    B = B[:-1]
    for i in range(n_terms):
        A[i,i] = -1
        if(i>0):
            A[i,i-1] = 1
    return np.linalg.inv(A)@B

# %%
def weno(f):
    eps = 1e-6
    length = len(f)
    coeff = [flux_coeff(1,2),flux_coeff(0,2)]
    gamma = np.array([1/3,2/3])

    fLplus = np.dot(f[1:3],coeff[0])
    fRplus = np.dot(f[2:4],coeff[1])
    fLminus= np.dot(f[0:2],coeff[0])
    fRminus = np.dot(f[1:3],coeff[1])

    beta_plus = np.array([f[2]-f[1],f[3]-f[2]])**2
    alpha_plus = gamma / (eps + beta_plus)**2
    w_plus = alpha_plus / np.sum(alpha_plus)
    fPhalf = w_plus[0]*fLplus + w_plus[1]*fRplus

    beta_minus = np.array([f[1]-f[0],f[2]-f[1]])**2
    alpha_minus = gamma / (eps + beta_minus)**2
    w_minus = alpha_minus / np.sum(alpha_minus)
    fMhalf = w_minus[0]*fLminus + w_minus[1]*fRminus

    return fMhalf, fPhalf


# %%
def lax_fried(fs,u,dx,nx):

    fMhalf, fPhalf = weno(fs)
    uMhalf, uPhalf = weno(u)

    fMplus = 0.5*(fMhalf + np.abs(uMhalf)*u)
    fMminus = 0.5*(fMhalf - np.abs(uMhalf)*u)

    fPplus = 0.5*(fPhalf + np.abs(uPhalf)*u)
    fPminus = 0.5*(fPhalf - np.abs(uPhalf)*u)

    zeros = np.zeros((nx))
    G = np.amin([zeros,u/np.abs(u)],axis=0)*(fPplus-fMplus)/dx\
           + np.amax([zeros,u/np.abs(u)],axis=0)*(fPminus-fMminus)/dx
    return -1 * G

# %%
def rungeKuta(u, dt, dx, nx):
    def f(u):
        return -1*u*u
    k1 = lax_fried(f(u),u,dx,nx)
    u1 = u + k1*dt/2
    k2 = lax_fried(f(u1),u1,dx,nx)
    u2 = u + k2*dt
    return u2

# %%
def main():
    dt = 0.1; tt = 0.9; nt = int(tt/dt)+1; t = np.linspace(0,tt,nt)
    l = 1; nx = 10+1; dx = l/(nx-1); x = np.linspace(0,l,nx)
    u = np.sin(2*np.pi*x) + 0.01

    ans = np.zeros((nt,nx))
    for it in range(nt-1):
        ans[it,:] = u
        u = rungeKuta(u,dt,dx,nx)
    ans[-1,:] = u
    print(ans)
    animate(x,ans)

    '''
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.interactive(False)
    plt.plot(x,u)
    plt.ylabel('some numbers')
    plt.ioff()
    plt.show()
    #import time
    #time.sleep(2)
    '''

# %%
if __name__ == "__main__":
    main()
