#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 01:38:50 2021

@author: Volodymyr
"""
from math import sqrt, pi,sin,cos
from scipy.constants import e, m_p, c, epsilon_0
from scipy.constants import physical_constants
import numpy as np
import copy
import random

class IBS_Martini_rates(object):

    def __init__(self, twiss,beam_dict,bunch,nu_x,nu_y,nu_z,log_c,fact):
        
        
        self.beam_dict=beam_dict
        self.twiss=twiss
        
        self.nu_x=nu_x
        self.nv_y=nu_y
        self.nu_z=nu_z
        
        self.log_c=log_c
    
    
        #Ion beam parameters

        self.Z = self.beam_dict['Z']
        self.mass =self.beam_dict['mass']
        self.KE=self.beam_dict['KE']

        self.k_ke = 8.9875517873681764E9; #Coulomb constant
        self.emit_x=self.beam_dict['emit_x']
        self.emit_y=self.beam_dict['emit_y']
        self.dp_p=self.beam_dict['dp_p']
        self.N_ptcl=self.beam_dict['N_ptcl']
        self.sigma_s0=self.beam_dict['sigma_s0']

        gamma = 1+self.KE/self.mass;
        beta = sqrt(gamma*gamma-1)/gamma;
        
        self.gamma=gamma
        self.beta=beta

        self.r_i = self.k_ke*self.Z**2*e*1e-6/self.mass; #Classical ion radius   
        self.bunched = bunch;

        self.energy_spread = beta*beta*self.dp_p;
        p0= gamma*self.mass*1e6*e*beta/c;
        
        ##Lattice functions
        s_v=self.twiss['s']
        betx_v=self.twiss['betx']
        alfx_v=self.twiss['alfx']
        mux_v=self.twiss['mux']
        dx_v=self.twiss['dxst']*fact
        dpx_v=self.twiss['dpxst']*fact
        bety_v=self.twiss['bety']
        alfy_v=self.twiss['alfy']
        muy_v=self.twiss['muy']
        dy_v=self.twiss['dyst']*fact
        dpy_v=self.twiss['dpyst']*fact
        
        l_element_v=[0]
        
        n_element=len(betx_v)
        
        for n in range(n_element-1):
            l_element_v.append(s_v[n+1]-s_v[n])
        
        circ=twiss.summary.length
        
        f0 = beta*c/circ; #revolution frequency.
        w0= 2*pi*f0; #angular frequency.
        
        ##Containers for 3D grid integration
        
        #Containers for 3D grid integration
        class TrigonometryStorageUV:
            sin_u2_cos_v2=g1=g2_1=g2_2=0
            
            def __init__(self,sin_u2_cos_v2,g1,g2_1,g2_2):
                self.sin_u2_cos_v2 = sin_u2_cos_v2
                self.g1 = g1
                self.g2_1 = g2_1
                self.g2_2 = g2_2 
                
        class TrigonometryStorageV:
            sin_v=cos_v=0
            
            def __init__(self,sin_v,cos_v ):
                self.sin_v = sin_v
                self.cos_v = cos_v
                
        class TrigonometryStorageU:
            sin_u=sin_u2=cos_u2=g3=0
            storage_uv=[]
            
            def __init__(self,sin_u,sin_u2,cos_u2,g3,storage_uv):
                self.sin_u = sin_u
                self.sin_u2 = sin_u2
                self.cos_u2 = cos_u2    
                self.g3 = g3
                self.storage_uv = storage_uv        
        
        class OpticalStorage:
            a=b2=c2=d2=dtld=k1=k2=k3=0
            
            def __init__(self,a,b2,c2,d2,dtld,k1,k2,k3):    
                self.a = a
                self.b2 = b2
                self.c2 = c2    
                self.d2 = d2
                self.dtld = dtld
                self.k1 = k1
                self.k2 = k2 
                self.k3 = k3
        
        rx_ibs=ry_ibs=rs_ibs=0
        
       ###Solver section###

        sigma_xbet=[None] * n_element
        sigma_y=[None] * n_element
        alfx2=[None] * n_element
        sigma_xp=[None] * n_element
        alfy2=[None] * n_element
        sigma_yp=[None] * n_element
        
        #Bunch properties at each element
        for i in range(n_element):
            sigma_xbet[i]=sqrt(betx_v[i]*self.emit_x);
            sigma_y[i]=sqrt(bety_v[i]*self.emit_y);
            alfx2[i]=alfx_v[i]**2;
            sigma_xp[i]=sqrt((1+alfx2[i])*self.emit_x/betx_v[i]);
            alfy2[i]=alfy_v[i]**2;
            sigma_yp[i]=sqrt((1+alfy2[i])*self.emit_y/bety_v[i]);
        
        #abcdk               
        os_all=[None] * n_element
        
        for i in range(n_element):
            betx = betx_v[i];
            alfx = alfx_v[i];
            dx = dx_v[i];
            dpx = dpx_v[i];
            alfy = alfy_v[i];
            d_tld = alfx*dx+betx*dpx;
            d_tld = alfx*dx+betx*dpx;
            sigma_x = sqrt(sigma_xbet[i]*sigma_xbet[i]+dx*dx*self.dp_p*self.dp_p);
            sigma_tmp = self.dp_p*sigma_xbet[i]/(gamma*sigma_x);
            q = 2*beta*gamma*sqrt(sigma_y[i]/self.r_i);
            a = sigma_tmp*sqrt(1+alfx*alfx)/sigma_xp[i];
            b2 = sigma_tmp*sqrt(1+alfy*alfy)/sigma_yp[i];
            b2 = b2**2;
            c2 = q*sigma_tmp;
            c2 = c2**2;
            d2 = self.dp_p*dx/sigma_x;
            d2 = d2**2;
            dtld = self.dp_p*d_tld/sigma_x;
            k1 = 1.0 / c2;
            k2 = a * a * k1;
            k3 = b2 * k1;
            os=OpticalStorage(a,b2,c2,d2,dtld,k1,k2,k3);
            os_all[i]=copy.deepcopy(os);
            os=None
        
        #Coef f

        storage_u=[]
        storage_v=[];
        dv = 2*pi/self.nv_y;
        v = -0.5*dv;
        
        for iv in range(self.nv_y):
            v += dv;
            storage_v.append(TrigonometryStorageV(sin(v), cos(v)));
            
        du = pi/nu_x;
        u = -0.5*du;
        for iu in range(nu_x):
            u += du;
            sin_u = sin(u);
            sin_u2 = sin_u **2;
            cos_u2 = 1 - sin_u2;
            g3 = 1 - 3 * cos_u2;
            storage_uv=[]
            for iv in range(self.nv_y):
                sin_u2_cos_v2 = sin_u2 * storage_v[iv].cos_v * storage_v[iv].cos_v;
                g1 = 1 - 3 * sin_u2_cos_v2;
                g2_1 = 1 - 3 * sin_u2 * storage_v[iv].sin_v * storage_v[iv].sin_v;
                g2_2 = 6 * sin_u * storage_v[iv].sin_v * storage_v[iv].cos_v;
                tempUV = TrigonometryStorageUV(sin_u2_cos_v2, g1, g2_1, g2_2);
                storage_uv.append(tempUV);
            storage_u.append(TrigonometryStorageU(sin_u, sin_u2, cos_u2, g3, storage_uv));

        
        f1=[None]*n_element
        f2=[None]*n_element
        f3=[None]*n_element
        
        if (log_c > 0):
            for ie in range(n_element):
                duvTimes2Logc = 2*pi*pi/(self.nu_x*self.nv_y) * 2 * self.log_c;
                tempf1=tempf2=tempf3=0;
                for iu in range(nu_x):
                    tu = storage_u[iu];
                    sum1=sum2=sum3=0;
                    for iv in range(self.nv_y):
                        tv = storage_v[iv];
                        tuv = tu.storage_uv[iv];
                        tmp = os_all[ie].a * tv.sin_v - os_all[ie].dtld * tv.cos_v;
                        tmp = tmp**2;
                        inv_int_z = (tuv.sin_u2_cos_v2 + tu.sin_u2 * tmp + os_all[ie].b2 * tu.cos_u2) * os_all[ie].k1;
                        sum1 += tuv.g1 / inv_int_z;
                        sum2 += (tuv.g2_1 + tuv.g2_2 * os_all[ie].dtld / os_all[ie].a) / inv_int_z;
                        sum3 += tu.g3 / inv_int_z;
        
                    tempf1 += tu.sin_u * sum1;
                    tempf2 += tu.sin_u * sum2;
                    tempf3 += tu.sin_u * sum3;
        
                f1[ie]=tempf1 * os_all[ie].k1 * duvTimes2Logc;
                f2[ie]=tempf2 * os_all[ie].k2 * duvTimes2Logc;
                f3[ie]=tempf3 * os_all[ie].k3 * duvTimes2Logc;
        
        # if slicing by nz, longer but almost the same
        
        # f1=[None]*n_element
        # f2=[None]*n_element
        # f3=[None]*n_element
        # for ie in range(10):    
            
        #     duv = 2*pi*pi/(nu_x*nv_y);
        #     tempf1 = tempf2 = tempf3 = 0;  
        #     for iu in range(nu_x):
        #         tu = storage_u[iu];
        #         sum1 = sum2 = sum3 = 0;
        #         for iv in range(nv_y):
        #             tv = storage_v[iv];
        #             tuv = tu.storage_uv[iv];
        
        #             tmp = os_all[ie].a * tv.sin_v - os_all[ie].dtld * tv.cos_v;
        #             tmp = tmp**2;
        #             d_uv = (tuv.sin_u2_cos_v2 + tu.sin_u2 * tmp + os_all[ie].b2 * tu.cos_u2) * os_all[ie].k1;
        #             int_z = 0;
        #             dz = 20/(d_uv*nu_z);
        #             z = -0.5*dz;
        #             for iz in range(nu_z):
        #                    z += dz;
        #                    int_z += np.exp(-d_uv*z)*np.log(1+z*z)*dz;
        
        #             sum1 += int_z * tuv.g1;
        #             sum2 += int_z * (tuv.g2_1 + tuv.g2_2 * os_all[ie].dtld / os_all[ie].a);
        #             sum3 += int_z * tu.g3;
        
        #         tempf1 += tu.sin_u * sum1;
        #         tempf2 += tu.sin_u * sum2;
        #         tempf3 += tu.sin_u * sum3;
        
        #     f1[ie] = tempf1 * os_all[ie].k1* duv;
        #     f2[ie] = tempf2 * os_all[ie].k2* duv;
        #     f3[ie] = tempf3 * os_all[ie].k3* duv;
        # print(f1[1])
        # print(f2[1])
        # print(f3[1])

        #koefficient a 
        liambda = self.N_ptcl/(2*sqrt(pi)*self.sigma_s0);
        beta3 = beta**3;
        gamma4 = gamma**4;
        koef_a =c*self.r_i**2*liambda/(16*pi*sqrt(pi)*self.dp_p*beta3*gamma4)/(self.emit_x*self.emit_y);  
        
        rsi=[None]*n_element
        rxi=[None]*n_element
        ryi=[None]*n_element
        #Rates by element
        a = koef_a
        rx = 0;
        ry = 0;
        rs = 0;
        n=1;#n=2 for coasting                
        inv_circ = 1/circ;
        
        for ie in range(n_element-1):
            l_element = l_element_v[ie+1];
            rsi[ie] = n*a*(1-os_all[ie].d2)*f1[ie]*l_element*inv_circ;
            rxi[ie] = a*(f2[ie]+f1[ie]*(os_all[ie].d2+os_all[ie].dtld*os_all[ie].dtld))*l_element*inv_circ;
            ryi[ie] = a*f3[ie]*l_element*inv_circ;
            rx += rxi[ie];
            ry += ryi[ie];
            rs += rsi[ie];
        
        
        rates=[rx,ry,rs]
        rates_fix=[]
        
        for r in rates:
            if np.abs([r])<1e-8:
                rates_fix.append(0)
            else:
                rates_fix.append(r)
               
        self.rsx=rates_fix[0]
        self.rsy=rates_fix[1]
        self.rss=rates_fix[2]
             
class IBS_BM_rates(object):

    def __init__(self, twiss,beam_dict,bunch,nu_x,nu_y,nu_z,log_c,fact):
        
        self.beam_dict=beam_dict
        self.twiss=twiss
        
        self.nu_x=nu_x
        self.nv_y=nu_y
        self.nu_z=nu_z
        
        self.log_c=log_c
    
    
        #Ion beam parameters

        self.Z = self.beam_dict['Z']
        self.mass =self.beam_dict['mass']
        self.KE=self.beam_dict['KE']

        self.k_ke = 8.9875517873681764E9; #Coulomb constant
        self.emit_x=self.beam_dict['emit_x']
        self.emit_y=self.beam_dict['emit_y']
        self.dp_p=self.beam_dict['dp_p']
        self.N_ptcl=self.beam_dict['N_ptcl']
        self.sigma_s0=self.beam_dict['sigma_s0']
        
        gamma = 1+self.KE/self.mass;
        beta = sqrt(gamma*gamma-1)/gamma;
        dp_p2= self.dp_p**2
        inv_dp_p2=1/dp_p2
        gamma2=gamma**2
        
        self.gamma=gamma
        self.beta=beta

        self.r_i = self.k_ke*self.Z**2*e*1e-6/self.mass; #Classical ion radius   
        self.bunched = bunch;

        self.energy_spread = beta*beta*self.dp_p;
        p0= gamma*self.mass*1e6*e*beta/c;
        
        ##Lattice functions
        s_v=self.twiss['s']
        betx_v=self.twiss['betx']
        alfx_v=self.twiss['alfx']
        mux_v=self.twiss['mux']
        dx_v=self.twiss['dxst']*fact
        dpx_v=self.twiss['dpxst']*fact
        bety_v=self.twiss['bety']
        alfy_v=self.twiss['alfy']
        muy_v=self.twiss['muy']
        dy_v=self.twiss['dyst']*fact
        dpy_v=self.twiss['dpyst']*fact
        
        l_element_v=[0]
        
        n_element=len(betx_v)
        
        for n in range(n_element-1):
            l_element_v.append(s_v[n+1]-s_v[n])
        
        circ=twiss.summary.length
        
        f0 = beta*c/circ; #revolution frequency.
        w0= 2*pi*f0; #angular frequency.       
        
        
        
        def rd( x, y, z):
            lolim = 6.E-51
            uplim = 1.0E+48
        
            if ( \
                x < 0.0 or \
                y < 0.0 or \
                x + y < lolim or \
                z < lolim or \
                uplim < x or \
                uplim < y or \
                uplim < z ):
                print ( '' )
                print ( 'RD - Error!' )
                print ( '  Invalid input arguments.' )
                print ( '  X = %g' % ( x ) )
                print ( '  Y = %g' % ( y ) )
                print ( '  Z = %g' % ( z ) )
                print ( '' )
                ierr = 1
                value = 0.0
                return value
        
            ierr = 0
            xn = x
            yn = y
            zn = z
            sigma = 0.0
            power4 = 1.0
        
            while ( True ):
        
                mu = ( xn + yn + 3.0 * zn ) * 0.2
                xndev = ( mu - xn ) / mu
                yndev = ( mu - yn ) / mu
                zndev = ( mu - zn ) / mu
                epslon = max ( abs ( xndev ), max ( abs ( yndev ), abs ( zndev ) ) )
        
                if ( epslon < 0.01 ):
                    c1 = 3.0 / 14.0
                    c2 = 1.0 / 6.0
                    c3 = 9.0 / 22.0
                    c4 = 3.0 / 26.0
                    ea = xndev * yndev
                    eb = zndev * zndev
                    ec = ea - eb
                    ed = ea - 6.0 * eb
                    ef = ed + ec + ec
                    s1 = ed * ( - c1 + 0.25 * c3 * ed - 1.5 * c4 * zndev * ef )
                    s2 = zndev  * ( c2 * ef + zndev * ( - c3 * ec + zndev * c4 * ea ) )
                    value = 3.0 * sigma  + power4 * ( 1.0 + s1 + s2 ) / ( mu * np.sqrt ( mu ) )
                    return value
        
                xnroot = np.sqrt ( xn )
                ynroot = np.sqrt ( yn )
                znroot = np.sqrt ( zn )
                lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
                sigma = sigma + power4 / ( znroot * ( zn + lamda ) )
                power4 = power4 * 0.25
                xn = ( xn + lamda ) * 0.25
                yn = ( yn + lamda ) * 0.25
                zn = ( zn + lamda ) * 0.25

        #BM optical storage
        class Kernels :
            psi=sx=sp=sxp=inv_sigma=0
    
            def __init__(self,psi,sx,sp,sxp,inv_sigma):    
                self.psi=psi;
                self.sx=sx;
                self.sp=sp;
                self.sxp=sxp;
                self.inv_sigma=inv_sigma;
        class OpticalStorageBM : #variables only depends on the TWISS parameters and the energy.
            phi=dx2=dx_betax_phi_2=sqrt_betay=gamma_phi_2=0
    
            def __init__(self,phi,dx2,dx_betax_phi_2,sqrt_betay,gamma_phi_2):
                self.phi=phi;
                self.dx2=dx2; #D_x * D_x
                self.dx_betax_phi_2=dx_betax_phi_2; # D_x * D_x / (beta_x * beta_x) + phi * phi
                self.sqrt_betay=sqrt_betay; # sqrt(beta_y)
                self.gamma_phi_2=gamma_phi_2; # gamma * gamma * phi * phi

       ###Solver section###

        sigma_xbet=[None] * n_element
        sigma_y=[None] * n_element
        alfx2=[None] * n_element
        sigma_xp=[None] * n_element
        alfy2=[None] * n_element
        sigma_yp=[None] * n_element
        
        #Bunch properties at each element
        for i in range(n_element):
            sigma_xbet[i]=sqrt(betx_v[i]*self.emit_x);
            sigma_y[i]=sqrt(bety_v[i]*self.emit_y);
            alfx2[i]=alfx_v[i]**2;
            sigma_xp[i]=sqrt((1+alfx2[i])*self.emit_x/betx_v[i]);
            alfy2[i]=alfy_v[i]**2;
            sigma_yp[i]=sqrt((1+alfy2[i])*self.emit_y/bety_v[i]);

        optical_strage=[None] * n_element
        for i in range(n_element):
            
            dx_betax = dx_v[i]/betx_v[i];
            dx2 = dx_v[i] * dx_v[i];
            phi = dpx_v[i]+alfx_v[i]*dx_betax;
            dx_betax_phi_2 = dx_betax*dx_betax + phi*phi;
            sqrt_betay = sqrt(bety_v[i]);
            gamma_phi_2 = gamma*gamma*phi*phi;
            optical_strage[i]= OpticalStorageBM(phi,dx2,dx_betax_phi_2,sqrt_betay,gamma_phi_2);
            
        kernels=[None] * n_element
        
        for i in range(n_element):
            betx = betx_v[i];
            bety = bety_v[i];
            dx = dx_v[i];
            dpx = dpx_v[i];
            sigma_x = sqrt(sigma_xbet[i]*sigma_xbet[i]+dx*dx*self.dp_p*self.dp_p);
            sigma_y = sqrt(bety_v[i]*self.emit_y);
            inv_sigma = 1/(sigma_x*sigma_y);
        
            ax = betx/self.emit_x;
            lamda_1 = bety/self.emit_y; #lamda_1 = ay.
            asa = ax*optical_strage[i].dx_betax_phi_2 + inv_dp_p2;
            a1 = gamma2*asa;
            a2 = (ax-a1)/2;
            a1 = (ax+a1)/2;
        
            lamda_2 = sqrt(a2*a2+ax*ax*optical_strage[i].gamma_phi_2);
            lamda_3 = a1 - lamda_2;
            tmp1 = 3/lamda_2;
            lamda_2 = a1 + lamda_2;
        
            inv_lamda_1 = 1/lamda_1;
            inv_lamda_2 = 1/lamda_2;
            inv_lamda_3 = 1/lamda_3;
        
            r1 = rd(inv_lamda_2, inv_lamda_3, inv_lamda_1);
            r2 = rd(inv_lamda_3, inv_lamda_1, inv_lamda_2);
            r3 = 3*sqrt(lamda_1*lamda_2*lamda_3)-r1-r2;
        
            r1 *= inv_lamda_1*2;
            r2 *= inv_lamda_2;
            r3 *= inv_lamda_3;
        
            psi = -r1 + r2 + r3;
        
            sxp = tmp1*ax*optical_strage[i].gamma_phi_2*(r3-r2);
            tmp1 = tmp1*a2;
            tmp2 = 1 + tmp1;
            tmp1 = 1 - tmp1;
            sp = gamma2*(r1-r2*tmp1-r3*tmp2)/2;
            sx = (r1-r2*tmp2-r3*tmp1)/2;
            kernels[i] = Kernels(psi,sx,sp,sxp,inv_sigma);
            
        lambd = 1/(2*sqrt(pi)*self.sigma_s0);

        beta3 = beta*beta*beta;
        gamma5 = gamma*gamma*gamma*gamma*gamma;
        
        lc=log_c
        c_bm=lambd*self.N_ptcl*self.r_i**2*c/(circ*6*sqrt(pi)*beta3*gamma5);
        c_bm=c_bm*lc;
        
        rsi=[None]*n_element
        rxi=[None]*n_element
        ryi=[None]*n_element
        
        rx = 0;
        ry = 0;
        rs = 0;
        n=1;
        
        for ie in range(n_element-1):
            l_element = l_element_v[ie+1];
            rxi[ie] = (betx_v[ie]*kernels[ie].inv_sigma*(kernels[ie].sx+optical_strage[ie].dx_betax_phi_2*kernels[ie].sp+kernels[ie].sxp)*l_element)*c_bm/self.emit_x;
            ryi[ie] = (bety_v[ie]*kernels[ie].inv_sigma*kernels[ie].psi*l_element)*c_bm/self.emit_y;
            rsi[ie] = (kernels[ie].inv_sigma*kernels[ie].sp*l_element)*n*c_bm/(dp_p2);
            rx += rxi[ie];
            ry += ryi[ie];
            rs += rsi[ie];  
            
        rates=[rx,ry,rs]
        rates_fix=[]
        
        for r in rates:
            if np.abs([r])<1e-8:
                rates_fix.append(0)
            else:
                rates_fix.append(r)
               
        self.rsx=rates_fix[0]
        self.rsy=rates_fix[1]
        self.rss=rates_fix[2]
    
class IBS_kick(object):
    
    def __init__(self, growthrate_x, growthrate_y, growthrate_z, dt,Dx=0,Dy=0):
        self.rxs=growthrate_x
        self.rys=growthrate_y
        self.rss=growthrate_z
        self.dt=dt
        self.Dx=Dx
        self.Dy=Dy
    
    def track(self, beam):
        
        betagamma=beam.betagamma
        
        emit_x=beam.epsn_x()/betagamma
        emit_y=beam.epsn_y()/betagamma
        
        num_sample=beam.macroparticlenumber
        
        betaX=1.971964#beam.beta_Twiss_x()
        betaY=3.41751#2#beam.beta_Twiss_y()
        SigDp=beam.sigma_dp()
        beta_z=beam.sigma_z()/beam.sigma_dp()
        
        
        if self.rxs>0:
            beam.xp +=np.random.normal(0,np.sqrt(2*self.rxs*self.dt*emit_x/betaX),num_sample)
        else:
            beam.xp *=np.exp(self.rxs*self.dt)
            
        if self.rys>0:
            beam.yp +=np.random.normal(0,np.sqrt(2*self.rys*self.dt*emit_y/betaY),num_sample)
        else:
            beam.yp *=np.exp(self.rys*self.dt)
        
        # subtract dispersion
        beam.x -= self.Dx* beam.dp
        beam.y -= self.Dy* beam.dp
        #kick
        if self.rss>0:
            dp_delta=np.random.normal(0,np.sqrt(2*self.rss*self.dt*SigDp**2),num_sample)
            beam.dp +=dp_delta
        else:
            beam.dp *=np.exp(self.rss*self.dt)
            
            
        # re-add dispersion
        beam.x += self.Dx* beam.dp
        beam.y += self.Dy* beam.dp
        