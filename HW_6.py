#Marie Calapa A732 HW6
#11/26/13
#Golf Ball Forces Simulation
#
import numpy as np
import matplotlib.pyplot as plt
import math

a = [0,0,-9.81] #grav
lofta = 0.35 #20 deg in rad (20deg = 0.35)
v0_1=33.53 #75 mph in 33.528m/s,90mph=40.23m/s,105mph=46.9m/s,120mph=53.53m/s
v0_2=40.23
v0_3=46.9
v0_4=53.53
mass=0.0459 #kgrams
#calc spin for all based on loft angle and club speed
spin = (5/7)*v0_1*math.sin(lofta)/0.02135
spin2 = (5/7)*v0_2*math.sin(lofta)/0.02135
spin3 = (5/7)*v0_3*math.sin(lofta)/0.02135
spin4 = (5/7)*v0_4*math.sin(lofta)/0.02135
#calc launch angl using loft angle
launch=np.arctan(5./(7+math.tan(lofta)**2)*math.tan(lofta))
#functions to calc coefficients of drag and lift given a spin and a velocity
def c_d(v,sp):
    if v != 0:
        S=0.02135*sp/v
        c_d = 0.225+0.2*S
    elif v == 0:
        c_d = 0
    return c_d

def c_l(v,sp):
    if v != 0:
        S=0.02135*sp/v
        #print(S)
        c_l = 0
        #c_l = 0.54*S**0.4
        #print(c_l)
    elif v == 0:
        c_l = 0
    return c_l

#calc drag and lift acceleration
def drag(v,sp):
    F = -c_d(v,sp)*1204*(math.pi*0.02135**2)*v/(2*mass)
    return F

def lift(v,sp):
    F_l = -c_l(v,sp)*1204*(math.pi*0.02135**2)*v/(2*mass)
    return F_l

aa=-0.087 #attack angle rads (-5 = -0.087)
la = aa+launch
dt=0.01
t=np.arange(0,5,dt)

#75
x=np.zeros(t.shape[0])
y=np.zeros(t.shape[0])
z=np.zeros(t.shape[0])
vx=np.zeros(t.shape[0])
vx[0]=v0_1*math.cos(la)
vy=np.zeros(t.shape[0])
vy[0]=0
vz=np.zeros(t.shape[0])
vz[0]=v0_1*math.sin(la)
#90
x2=np.zeros(t.shape[0])
y2=np.zeros(t.shape[0])
z2=np.zeros(t.shape[0])
vx2=np.zeros(t.shape[0])
vx2[0]=v0_2*math.cos(la)
vy2=np.zeros(t.shape[0])
vy2[0]=0
vz2=np.zeros(t.shape[0])
vz2[0]=v0_2*math.sin(la)
#105
x3=np.zeros(t.shape[0])
y3=np.zeros(t.shape[0])
z3=np.zeros(t.shape[0])
vx3=np.zeros(t.shape[0])
vx3[0]=v0_3*math.cos(la)
vy3=np.zeros(t.shape[0])
vy3[0]=0
vz3=np.zeros(t.shape[0])
vz3[0]=v0_3*math.sin(la)
#120
x4=np.zeros(t.shape[0])
y4=np.zeros(t.shape[0])
z4=np.zeros(t.shape[0])
vx4=np.zeros(t.shape[0])
vx4[0]=v0_4*math.cos(la)
vy4=np.zeros(t.shape[0])
vy4[0]=0
vz4=np.zeros(t.shape[0])
vz4[0]=v0_4*math.sin(la)


#RK in all directions (p,q = x,vx n,o=y,vy l,m=z,vz)
def rk(p,q,n,o,l,m,t,sp):
    dt=0.01
    
    for i in range(0,t.shape[0]-1,1):
        kvx_1 = dt*(q[0]+drag(q[i],sp)*t[i]+lift(m[i],sp)*t[i])
        kvy_1 = dt*(o[0]+drag(o[i],sp)*t[i])
        kvz_1 = dt*(a[2]*(t[i])+drag(m[i],sp)*t[i]+lift(q[i],sp)*t[i])
        kvx_2 = dt*(0.5*kvx_1 +drag(q[i],sp)*(t[i]+dt/2.)+lift(m[i],sp)*(t[i]+dt/2.))
        kvy_2 = dt*(0.5*kvy_1 +drag(o[i],sp)*(t[i]+dt/2.))
        kvz_2 = dt*(0.5*kvz_1 + a[2]*(t[i]+dt/2.)+drag(m[i],sp)*(t[i]+dt/2.)+lift(q[i],sp)*(t[i]+dt/2.))
        kvx_3 = dt*(0.5*kvx_2+drag(q[i],sp)*(t[i]+dt/2.)+lift(m[i],sp)*(t[i]+dt/2.))
        kvy_3 = dt*(0.5*kvy_2+drag(o[i],sp)*(t[i]+dt/2.))
        kvz_3 = dt*(0.5*kvz_2 + a[2]*(t[i]+dt/2.)+drag(m[i],sp)*(t[i]+dt/2.)+lift(q[i],sp)*(t[i]+dt/2.))
        kvx_4 = dt*(kvx_3+drag(q[i],sp)*(t[i]+dt)+lift(m[i],sp)*(t[i]+dt))
        kvy_4 = dt*(kvy_3+drag(o[i],sp)*(t[i]+dt))
        kvz_4 = dt*(kvz_3 + a[2]*(t[i]+dt)+drag(m[i],sp)*(t[i]+dt)+lift(q[i],sp)*(t[i]+dt))
        q[i+1] = q[i]+(1./6.)*kvx_1 + (1./3.)*kvx_2 + (1./3.)*kvx_3 + (1./6.)*kvx_4
        o[i+1] = o[i]+(1./6.)*kvy_1 + (1./3.)*kvy_2 + (1./3.)*kvy_3 + (1./6.)*kvy_4
        m[i+1] = m[i]+(1./6.)*kvz_1 + (1./3.)*kvz_2 + (1./3.)*kvz_3 + (1./6.)*kvz_4
        kx_1 = dt*(q[i]*t[i])
        ky_1 = dt*(o[i]*t[i])
        kz_1 = dt*(m[i]*t[i])
        kx_2 = dt*(0.5*kx_1+q[i]*(t[i]+dt/2.))
        ky_2 = dt*(0.5*ky_1+o[i]*(t[i]+dt/2.))
        kz_2 = dt*(0.5*kz_1+m[i]*(t[i]+dt/2.))
        kx_3 = dt*(0.5*kx_2+q[i]*(t[i]+dt/2.))
        ky_3 = dt*(0.5*ky_2+o[i]*(t[i]+dt/2.))
        kz_3 = dt*(0.5*kz_2+m[i]*(t[i]+dt/2.))
        kx_4 = dt*(kx_3+q[i]*(t[i]+dt))
        ky_4 = dt*(ky_3+o[i]*(t[i]+dt))
        kz_4 = dt*(kz_3+m[i]*(t[i]+dt))
        p[i+1] = p[i]+(1./6.)*kx_1 + (1./3.)*kx_2 + (1./3.)*kx_3 + (1./6.)*kx_4
        n[i+1] = n[i]+(1./6.)*ky_1 + (1./3.)*ky_2 + (1./3.)*ky_3 + (1./6.)*ky_4
        l[i+1] = l[i]+(1./6.)*kz_1 + (1./3.)*kz_2 + (1./3.)*kz_3 + (1./6.)*kz_4
    
    return p,q,n,o,l,m

#call it for all club speeds
x_fin,vx_fin,y_fin,vy_fin,z_fin,vz_fin = rk(x,vx,y,vy,z,vz,t,spin)
x_fin2,vx_fin2,y_fin2,vy_fin2,z_fin2,vz_fin2 = rk(x2,vx2,y2,vy2,z2,vz2,t,spin2)
x_fin3,vx_fin3,y_fin3,vy_fin3,z_fin3,vz_fin3 = rk(x3,vx3,y3,vy3,z3,vz3,t,spin3)
x_fin4,vx_fin4,y_fin4,vy_fin4,z_fin4,vz_fin4 = rk(x4,vx4,y4,vy4,z4,vz4,t,spin4)


#plot results for all velocities
plt.plot(x_fin,z_fin,label='v0 = 75 mph')
plt.plot(x_fin2,z_fin2,label='v0 = 90 mph')
plt.plot(x_fin3,z_fin3,label='v0 = 105 mph')
plt.plot(x_fin4,z_fin4,label='v0 = 120 mph')
plt.legend()
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')
plt.ylim(0,5)
plt.xlim(0,20)
plt.show()

