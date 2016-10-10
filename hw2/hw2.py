from pylab import *

# Prints answers to problem 1
print 'Column density = 3.09e20 cm^-2.'
print '1(a) 3.24e-24 cm^2.'
print '1(b) 3.24e-21 cm^2.'
print '1(c) 3.24e-18 cm^2.'

pctocm = 3.086e+18 # pc to cm
D = 100. * pctocm # Depth D in cm
ds = 0.01 * pctocm # step size of D in cm
n = 1. # number density of medium in cm^-3
nstep = D/ds+1 # number of steps of D

def solve_nu_eq(sigma_nu, I_nu0, S_nu):
    # functino to solve radiative transfer equation
    # takes sigma_nu, I_nu0 and S_nu as input
    s = linspace(0,D,num=nstep) # creates array of distance s
    i_nu = zeros((nstep, size(I_nu0))) # creates array of specific intensity along s
    i_nu[0,:] = I_nu0 # sets initial specific intensity
    alpha_nu = n * sigma_nu # calculates absorption coefficient
    for i in range(len(s)-1):
        i_nu[i+1] = i_nu[i] + alpha_nu * (S_nu - i_nu[i]) * ds # solves radiatio transfer equation
    return i_nu
    
def gen_sigma_nu(nu, sigma_func, args):
    # function to generate sigma_nu from a function
    # takes array of frequencies nu, function representing sigma sigma_func, 
    # and arguments args passed to sigma_func as input
    return sigma_func(nu, *args)
    
def gaussian(nu,A,c,nu_0):
    # Gaussian function with normalization A, sigma c and centroid nu_0
    return A*exp(-(nu-nu_0)**2./(2.*c**2.))

def const(nu, c):
    # constant function
    return c * ones(len(nu))

nu = linspace(-5,5) # creates array of frequencies

figure(figsize=(6,10))

# (a)
I_nu0 = ones(len(nu)) # I_nu0 = 1
sigma_nu = gen_sigma_nu(nu, const, (3.24e-18,)) # sigma = 3.24e-18
S_nu = 0.5

I_nu = solve_nu_eq(sigma_nu, I_nu0, S_nu)
subplot(611)
plot(nu, I_nu[-1])

ylabel(r'$I_\nu$')

# (b) 
I_nu0 = zeros(len(nu)) # I_nu0 = 0
sigma_nu = gen_sigma_nu(nu, gaussian, (1e-21, 1., 0.))
S_nu = 0.5

I_nu = solve_nu_eq(sigma_nu, I_nu0, S_nu)
subplot(612)
plot(nu, I_nu[-1])

ylabel(r'$I_\nu$')

# (c)
I_nu0 = ones(len(nu))
sigma_nu = gen_sigma_nu(nu, gaussian, (1e-21, 1., 0.))
S_nu = 2.

I_nu = solve_nu_eq(sigma_nu, I_nu0, S_nu)
subplot(613)
plot(nu, I_nu[-1])

ylabel(r'$I_\nu$')

# (d)
I_nu0 = ones(len(nu))
sigma_nu = gen_sigma_nu(nu, gaussian, (1e-21, 1., 0.))
S_nu = 0.5

I_nu = solve_nu_eq(sigma_nu, I_nu0, S_nu)
subplot(614)
plot(nu, I_nu[-1])

ylabel(r'$I_\nu$')

# (e)
I_nu0 = ones(len(nu))
sigma_nu = gen_sigma_nu(nu, gaussian, (5e-20, 1., 0.))
S_nu = 2.

I_nu = solve_nu_eq(sigma_nu, I_nu0, S_nu)
subplot(615)
plot(nu, I_nu[-1])

ylabel(r'$I_\nu$')

# (f)
I_nu0 = ones(len(nu))
sigma_nu = gen_sigma_nu(nu, gaussian, (5e-20, 1., 0.))
S_nu = 0.5

I_nu = solve_nu_eq(sigma_nu, I_nu0, S_nu)
subplot(616)
plot(nu, I_nu[-1])

ylabel(r'$I_\nu$')
xlabel(r'$\nu-\nu_0$ (Arbitrary units)')

tight_layout()

show()