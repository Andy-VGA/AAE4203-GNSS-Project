#!/usr/bin/env python3
"""Independent GPS L1 broadcast-ephemeris WLS-SPP benchmark baseline.

Requires numpy; does NOT use RTKLIB or MatRTKLIB. Uses GPS C1C pseudorange,
GPS broadcast ephemerides, Earth-rotation correction, satellite clock including
relativity/TGD, Klobuchar ionosphere and Saastamoinen-style troposphere.
Outputs only independently solved epochs (no interpolation or invented fixes).
This is an experimental baseline, not a verified high-accuracy solution.
"""
import argparse, csv, math
from datetime import datetime
from pathlib import Path
from collections import Counter
import numpy as np

C=299792458.0
MU=3.986005e14
OMEGA=7.2921151467e-5
F=-4.442807633e-10
GPS0=datetime(1980,1,6)
A=6378137.0
E2=6.69437999014e-3

def gps_time(year,mo,day,hh,mm,ss):
    return (datetime(year,mo,day,hh,mm)-GPS0).total_seconds()+ss

def fnum(s):
    try: return float(s.replace('D','E').replace('d','e'))
    except (ValueError,TypeError): return float('nan')

def parse_nav(path):
    with open(path,encoding='ascii',errors='replace') as f: lines=f.read().splitlines()
    header_end=next(i for i,l in enumerate(lines) if 'END OF HEADER' in l)
    alpha=[0.,0.,0.,0.]; beta=[0.,0.,0.,0.]
    for line in lines[:header_end]:
        if 'IONOSPHERIC CORR' in line and line[:4] in ('GPSA','GPSB'):
            vals=[fnum(line[5+i*12:5+(i+1)*12]) for i in range(4)]
            if line[:4]=='GPSA': alpha=vals
            else: beta=vals
    nav={}
    i=header_end+1
    while i<len(lines):
        line=lines[i]
        if len(line)<3 or not line[0].isalpha(): i+=1; continue
        system=line[0]
        count=4 if system in ('R','S') else 8
        block=lines[i:i+count]
        i+=count
        if system!='G' or len(block)!=8: continue
        try:
            sat=line[:3].strip()
            fields=line[4:23].split()
            yy,mo,day,hh,mm=map(int,fields[:5]); sec=float(fields[5])
            toc=gps_time(yy,mo,day,hh,mm,sec)
            clock=[fnum(line[23+k*19:23+(k+1)*19]) for k in range(3)]
            vals=[]
            for b in block[1:]:
                vals.extend(fnum(b[4+k*19:4+(k+1)*19]) for k in range(4))
            if not np.isfinite(vals[7]) or not np.isfinite(vals[8]): continue
            week=int(vals[18]) if np.isfinite(vals[18]) else int(toc//604800)
            toe=week*604800+vals[8]
            nav.setdefault(sat,[]).append(dict(toc=toc,toe=toe,af0=clock[0],af1=clock[1],af2=clock[2],v=vals))
        except Exception: continue
    for records in nav.values(): records.sort(key=lambda x:x['toe'])
    return nav,alpha,beta

def parse_obs(path):
    with open(path,encoding='ascii',errors='replace') as f: lines=f.read().splitlines()
    head=next(i for i,l in enumerate(lines) if 'END OF HEADER' in l)
    types={}
    for line in lines[:head]:
        if 'SYS / # / OBS TYPES' in line and line[0].isalpha():
            system=line[0]; n=int(line[3:6]); types[system]=line[7:60].split()[:n]
    gps_types=types['G']; ci=gps_types.index('C1C'); si=gps_types.index('S1C') if 'S1C' in gps_types else None
    i=head+1
    while i<len(lines):
        l=lines[i];i+=1
        if not l.startswith('>'): continue
        try:
            fields=l[1:29].split()
            year,mo,day,hh,mm=map(int,fields[:5]); sec=float(fields[5]); n=int(l[32:35]); flag=int(l[30:32])
            t=gps_time(year,mo,day,hh,mm,sec)
        except (ValueError,IndexError): continue
        obs=[]
        for _ in range(n):
            if i>=len(lines): break
            s=lines[i];i+=1
            if s[:1]!='G': continue
            sid=f'G{int(s[1:3]):02d}'; fields=[s[3+k*16:3+(k+1)*16] for k in range(len(gps_types))]
            pr=fnum(fields[ci][:14]) if ci<len(fields) else float('nan')
            snr=fnum(fields[si][:14]) if si is not None and si<len(fields) else float('nan')
            if np.isfinite(pr) and 1e7<pr<5e7 and (not np.isfinite(snr) or snr>=15):
                obs.append((sid,pr,snr))
        if flag in (0,1):yield t,obs

def wrap(t):
    if t>302400: t-=604800
    elif t<-302400: t+=604800
    return t

def select_eph(nav,sat,t):
    if sat not in nav:return None
    return min(nav[sat],key=lambda x:abs(wrap(t-x['toe'])))

def satellite(eph,t):
    v=eph['v']
    crs,dn,m0,cuc,e,cus,sqrtA,toe,cic,omega0,cis,i0,crc,omega,omega_dot,idot=(
        v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16])
    a=sqrtA*sqrtA
    tk=wrap(t-eph['toe'])
    n0=math.sqrt(MU/a**3)
    m=m0+(n0+dn)*tk
    E=m
    for _ in range(15):
        En=m+e*math.sin(E)
        if abs(En-E)<1e-13:break
        E=En
    sinE=math.sin(E); cosE=math.cos(E)
    nu=math.atan2(math.sqrt(1-e*e)*sinE,cosE-e)
    phi=nu+omega
    u=phi+cuc*math.cos(2*phi)+cus*math.sin(2*phi)
    r=a*(1-e*cosE)+crc*math.cos(2*phi)+crs*math.sin(2*phi)
    inc=i0+idot*tk+cic*math.cos(2*phi)+cis*math.sin(2*phi)
    xp=r*math.cos(u); yp=r*math.sin(u)
    bigO=omega0+(omega_dot-OMEGA)*tk-OMEGA*toe
    cb=math.cos(bigO);sb=math.sin(bigO);ci=math.cos(inc);si=math.sin(inc)
    xyz=np.array([xp*cb-yp*ci*sb,xp*sb+yp*ci*cb,yp*si])
    dt=t-eph['toc']
    clock=eph['af0']+eph['af1']*dt+eph['af2']*dt*dt+F*e*sqrtA*sinE-v[22]
    return xyz,clock

def llh(xyz):
    x,y,z=xyz
    lon=math.atan2(y,x); p=math.hypot(x,y)
    lat=math.atan2(z,p*(1-E2))
    for _ in range(5):
        N=A/math.sqrt(1-E2*math.sin(lat)**2)
        h=p/math.cos(lat)-N
        lat=math.atan2(z,p*(1-E2*N/(N+h)))
    N=A/math.sqrt(1-E2*math.sin(lat)**2)
    h=p/math.cos(lat)-N
    return lat,lon,h

def enu_matrix(lat,lon):
    sl,cl=math.sin(lat),math.cos(lat);so,co=math.sin(lon),math.cos(lon)
    return np.array([[-so,co,0],[-sl*co,-sl*so,cl],[cl*co,cl*so,sl]])

def iono_delay(t,lat,lon,az,el,alpha,beta):
    # GPS ICD Klobuchar: semicircles; delay in metres
    phi=lat/math.pi; lam=lon/math.pi
    E=el/math.pi
    psi=0.0137/(E+0.11)-0.022
    phii=max(-0.416,min(0.416,phi+psi*math.cos(az)))
    lami=lam+psi*math.sin(az)/math.cos(phii*math.pi)
    phim=phii+0.064*math.cos((lami-1.617)*math.pi)
    local=(43200*lami+t%86400)%86400
    amp=max(0,sum(alpha[i]*phim**i for i in range(4)))
    per=max(72000,sum(beta[i]*phim**i for i in range(4)))
    xx=2*math.pi*(local-50400)/per
    factor=1+16*(0.53-E)**3
    return C*factor*(5e-9+(amp*(1-xx*xx/2+xx**4/24) if abs(xx)<1.57 else 0))

def tropo_delay(el,h):
    # Standard atmosphere + simple zenith hydrostatic/wet delay model
    height=max(-100,min(10000,h))
    pres=1013.25*(1-2.2557e-5*height)**5.2568
    temp=15.0-0.0065*height+273.15
    vapor=6.108*0.7*math.exp((17.15*(temp-273.15)-4684)/(temp-38.45)) if temp>200 else 0
    zhd=0.0022768*pres
    zwd=0.002277*(1255/temp+0.05)*vapor
    return (zhd+zwd)/max(0.1,math.sin(el))

def build_measurements(t,obs,nav,xyz,alpha,beta,mask=15):
    lat,lon,h=llh(xyz)
    rot=enu_matrix(lat,lon)
    arr=[]
    for sid,pr,snr in obs:
        eph=select_eph(nav,sid,t)
        if eph is None or abs(wrap(t-eph['toe']))>4*3600:continue
        if np.isfinite(eph['v'][21]) and eph['v'][21]!=0:continue
        ts=t-pr/C
        for _ in range(2):
            sat,clk=satellite(eph,ts)
            ts=t-pr/C+clk
        sat,clk=satellite(eph,ts)
        tau=pr/C
        cos=math.cos(OMEGA*tau); sin=math.sin(OMEGA*tau)
        # Earth rotation compensation during propagation
        sat=np.array([cos*sat[0]+sin*sat[1],-sin*sat[0]+cos*sat[1],sat[2]])
        dr=sat-xyz
        rho=np.linalg.norm(dr)
        local=rot@dr
        el=math.atan2(local[2],math.hypot(local[0],local[1]))
        if math.degrees(el)<mask:continue
        az=math.atan2(local[0],local[1])%(2*math.pi)
        ion=iono_delay(t,lat,lon,az,el,alpha,beta)
        trop=tropo_delay(el,h)
        corrected=pr+C*clk-ion-trop
        weight=(max(0.1,math.sin(el))**2)* (min(1.0,max(0.3,snr/40))**2 if np.isfinite(snr) else 1)
        arr.append((sid,sat,rho,corrected,weight,math.degrees(el)))
    return arr

def solve_epoch(t,obs,nav,init,alpha,beta,mask=15):
    xyz=init[:3].copy(); bias=float(init[3]); last=None
    for iteration in range(12):
        measurements=build_measurements(t,obs,nav,xyz,alpha,beta,mask)
        if len(measurements)<4:return None
        H=[];v=[];w=[]
        for sid,sat,rho,pr,wt,el in measurements:
            H.append(np.r_[-(sat-xyz)/rho,1.])
            v.append(pr-rho-bias)
            w.append(wt)
        H=np.array(H);v=np.array(v);w=np.array(w)
        try:
            dw=np.linalg.lstsq(H*np.sqrt(w[:,None]),v*np.sqrt(w),rcond=None)[0]
        except np.linalg.LinAlgError:return None
        if not np.isfinite(dw).all() or np.linalg.norm(dw[:3])>2e5:return None
        xyz+=dw[:3];bias+=dw[3]
        if np.linalg.norm(dw[:3])<0.001:
            res=v-H@dw
            return dict(xyz=xyz,bias=bias,ns=len(H),iterations=iteration+1,rms=float(np.sqrt(np.mean(res**2))))
    return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--obs',default='urbannav_deep1.obs')
    ap.add_argument('--nav',default='urbannav_deep1.nav')
    ap.add_argument('--out',default='GXX_benchmark.csv')
    ap.add_argument('--mask',type=float,default=15.)
    args=ap.parse_args()
    nav,alpha,beta=parse_nav(args.nav)
    print('GPS broadcast nav satellites:',len(nav),'GPSA',alpha,'GPSB',beta,flush=True)
    rows=[];allrows=[];status=Counter()
    init=np.array([-2418000.,5385000.,2406000.,0.])
    for k,(t,obs) in enumerate(parse_obs(args.obs)):
        sol=solve_epoch(t,obs,nav,init,alpha,beta,args.mask)
        week=int(t//604800);tow=t-week*604800
        if sol is None:
            status['failed']+=1
            allrows.append([week,f'{tow:.3f}','nan','nan','nan',len(obs),'failed'])
        else:
            xyz=sol['xyz']
            init=np.r_[xyz,sol['bias']]
            rows.append([week,f'{tow:.3f}',*[f'{p:.3f}' for p in xyz]])
            allrows.append([week,f'{tow:.3f}',*[f'{p:.3f}' for p in xyz],sol['ns'],'solved'])
            status['solved']+=1
        if (k+1)%300==0:print('epochs',k+1,dict(status),flush=True)
    with open(args.out,'w',newline='') as f:
        w=csv.writer(f);w.writerow(['gps_week','gps_tow_s','x_ecef_m','y_ecef_m','z_ecef_m']);w.writerows(rows)
    diag=Path(args.out).with_name('urbannav_benchmark_diagnostic.csv')
    with open(diag,'w',newline='') as f:
        w=csv.writer(f);w.writerow(['gps_week','gps_tow_s','x_ecef_m','y_ecef_m','z_ecef_m','num_satellites','status']);w.writerows(allrows)
    print('Total epochs:',len(allrows),'solved:',len(rows),'availability:',len(rows)/len(allrows)*100,'%',flush=True)
    print('Submission candidate:',args.out,'diagnostics:',diag,flush=True)

if __name__=='__main__':main()
