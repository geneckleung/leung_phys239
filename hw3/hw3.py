from pylab import *

e = 1.60217662e-19
Z = 1.
me = 9.10938356e-31
a0 = 5.29e-11
epsilon = 8.85418782e-12
c = 3e8

def solve(x0,y0,vx0,vy0,n,dt):
    # calculates trajectory
    # takes x0, y0 in units of Bohr radius, vx, vy in units of cm/s, number of steps n and step size dt
    # create arrays
    x = zeros(n+1)
    y = zeros(n+1)
    vx = zeros(n+1)
    vy = zeros(n+1)
    ax = zeros(n+1)
    ay = zeros(n+1)
    t = zeros(n+1)
    T = zeros(n+1)
    U = zeros(n+1)
    E = zeros(n+1)
    # set initial conditions
    x[0] = x0 * a0
    y[0] = y0 * a0
    vx[0] = vx0 / 1e2
    vy[0] = vy0 / 1e2
    ax[0] = -Z*e**2.*x[0]/(x[0]**2.+y[0]**2.)**(3./2.)/me/4./pi/epsilon
    ay[0] = -Z*e**2.*y[0]/(x[0]**2.+y[0]**2.)**(3./2.)/me/4./pi/epsilon
    T[0] = 1./2.*me*(vx[0]**2.+vy[0]**2.)# + 3./8.*me*(vx[0]**2.+vy[0]**2.)**2./c**2.
    U[0] = -Z*e**2./(x[0]**2.+y[0]**2.)**(1./2.)/4./pi/epsilon
    E[0] = T[0] + U[0]
    print E[0]
    # main loop
    for i in range(n):
        # breaks loop if outside region of interest
        if abs(x[i])>1000.*a0 or abs(y[i])>1000.*a0:
            print 'Number of steps = ', i
            end = i+1
            break
        t[i+1] = t[i] + dt
        vx[i+1] = vx[i] + ax[i]*dt
        vy[i+1] = vy[i] + ay[i]*dt
        x[i+1] = x[i] + vx[i]*dt
        y[i+1] = y[i] + vy[i]*dt
        ax[i+1] = -Z*e**2.*x[i+1]/(x[i+1]**2.+y[i+1]**2.)**(3./2.)/me/4./pi/epsilon
        ay[i+1] = -Z*e**2.*y[i+1]/(x[i+1]**2.+y[i+1]**2.)**(3./2.)/me/4./pi/epsilon
        T[i+1] = 1./2.*me*(vx[i+1]**2.+vy[i+1]**2.)# + 3./8.*me*(vx[i+1]**2.+vy[i+1]**2.)**2./c**2.
        U[i+1] = -Z*e**2./(x[i+1]**2.+y[i+1]**2.)**(1./2.)/4./pi/epsilon
        E[i+1] = T[i+1] + U[i+1]
        end = i+1
    # trim array if necessary
    x = x[:end] / a0
    y = y[:end] / a0
    vx = vx[:end] * 1e2
    vy = vy[:end] * 1e2
    ax = ax[:end]
    ay = ay[:end]
    t = t[:end]
    E = E[:end]
    U = U[:end]
    T = T[:end]
    return (x,y,vx,vy,ax,ay,t,E,U,T)

    
def analyze(x0,y0,vx0,vy0,n,dt,draw=False):
    # makes plots of trajectory, v and a; performs Fourier transform to get spectrum, returns peak frequency
    # takes same arguments as solve
    x,y,vx,vy,ax,ay,t,E,U,T=solve(x0,y0,vx0,vy0,n,dt)
    
    print 'Max. t = ', max(t)
    print 'Max. v/c = ', max(sqrt(vx**2+vy**2))/3e10
    
    # performs Fourier transform of acceleration
    frq = fftfreq(t.shape[-1],d=dt)
    fx = fft(ax)
    fy = fft(ay)
    spec = absolute(fx)**2.+absolute(fy)**2.
    peak = abs(frq[argmax(spec)])
    print 'Peak = ', peak
    
    # make plots
    if draw:
        # plot trajectory
        figure(figsize=(6,6))
        plot(x,y)
        plot(0,0,'ko')
        xlim(-1000,1000)
        ylim(-1000,1000)
        xlabel('x (a_0)')
        ylabel('y (a_0)')
        caption = ' x_0 = '+str(x0)+'\n y_0 = '+str(y0)+'\n v_0 = %.1e' % vx0
        figtext(0.65,0.8,caption)
        
        # plot v and a
        figure()
        subplot(221)
        plot(t,vx)
        xlabel('t (s)')
        ylabel('v_x (cm/s)')
        subplot(222)
        plot(t,vy)
        xlabel('t (s)')
        ylabel('v_y (cm/s)')
        subplot(223)
        plot(t,ax)
        xlabel('t (s)')
        ylabel('a_x (cm/s^2)')
        subplot(224)
        plot(t,ay)
        xlabel('t (s)')
        ylabel('a_y (cm/s^2)')
        tight_layout()

        # plot Fourier transform
        figure()
        ax1 = subplot2grid((2,2),(0,0))
        ax2 = subplot2grid((2,2),(0,1))
        ax3 = subplot2grid((2,2),(1,0),colspan=2)
        ax1.plot(frq,fx.real,label='Real part')
        ax1.plot(frq,fx.imag,label='Imaginary part')
        ax1.plot(frq,absolute(fx),label='Absolute value')
        ax1.legend(loc=4,fontsize=8)
        ax1.set_xlim(0,0.005*max(frq))
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Power')
        ax2.plot(frq,fy.real,label='Real part')
        ax2.plot(frq,fy.imag,label='Imaginary part')
        ax2.plot(frq,absolute(fy),label='Absolute value')
        ax2.legend(loc=4,fontsize=8)
        ax2.set_xlim(0,0.005*max(frq))
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('Power')
        ax3.plot(frq,spec)
        ax3.set_xlim(0,0.005*max(frq))
        ax3.set_xlabel('Frequency (Hz)')
        ax3.set_ylabel('Power')
        ax3.axvline(x=peak,color='k',ls='--')
        tight_layout()
        figtext(0.35,0.9,'FFT of a_x')
        figtext(0.85,0.9,'FFT of a_y')
        figtext(0.6,0.4,'Squared sum of FFT of a_x and a_y')
    return peak
    
p1 = analyze(-999., 80., 1e8, 0., int(1e5),5e-18)
p2 = analyze(-999., 60., 1e8, 0., int(1e5),5e-18,draw=True)
p3 = analyze(-999., 40., 1e8, 0., int(1e5),5e-18)
p4 = analyze(-999., 20., 1e8, 0., int(1e5),5e-18)


p5 = analyze(-999., 50., 3e8, 0., int(1e5),5e-18)
p6 = analyze(-999., 50., 2e8, 0., int(1e5),5e-18)
p7 = analyze(-999., 50., 1e8, 0., int(1e5),5e-18)
p8 = analyze(-999., 50., 5e7, 0., int(1e5),5e-18)

figure()
plot([80,60,40,20],[p1,p2,p3,p4])
xlabel('Impact parameter (a_0)')
ylabel('Peak frequency (Hz)')
figtext(0.7,0.85,'v_0 = 1e8 cm/s')

figure()
plot([3e8,2e8,1e8,5e7],[p5,p6,p7,p8])
xlabel('Initial velocity (cm/s)')
ylabel('Peak frequency (Hz)')
figtext(0.7,0.85,'b = 50*a_0')
show()
