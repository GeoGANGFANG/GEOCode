/* 2-D Lowrank Finite-difference wave extrapolation */
/*
  Copyright (C) 2015 China Geology Survey
  Author: Gang Fang
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <rsf.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "source.h"

static float ***G;
static int *sx, *sz;
static int LEN;
static marg;

static float fd(float **data, int ix, int iz)
/*< Lowrank finite difference >*/
{
    float res = 0.0;
    int il;
    for (il=0; il<size; il++) {
	res += 0.5*( data[ix-sx[il]][iz-sz[il]] + data[ix+sx[il]][iz+sz[il]] )*G[il][ix][iz];
    }
    return res;
}

int main(int argc, char* argv[])
{
    clock_t tstart, tend;
    double duration;
    
    /* omp */
    int nth;
    
    /* flag */
    bool verb;
    
    /* I/O */
    sf_file fvel, fsource, fwf, frec;
    sf_file fG, fsx, fsz;
    
    sf_axis at, ax, az;
    
    /* I/O arrays */
    float *src, **record;
    float **vel;
    
    float *sxtmp, *sztmp;
    
    /* grid indec variables */
    int nx, nz, nxz, nt, ix, iz, it;
    int nxb, nzb;
    float dt, dx, dz;
    float ox, oz;
    
    /* caculate arrays */
    float **pn1, **pn0;

    /* source */
    spara sp={0};
    bool  srcdecay;
    int   srcrange;
    float srctrunc;
    float slx, slz;
    int   spx, spz;
    
    /* PML */
    int mx, mz; /* margin */

    /* options */
    bool freesurface;
    float gdep;
    int gp;
    int snapinter;

    tstart = clock();
    sf_init(argc, argv);
    if (!sf_getbool("verb", &verb)) verb=false; /*verbosity*/
    
    /* set I/O file */
    fsource = sf_input("in");  /*source wavelet*/
    fvel    = sf_input("vel"); /*velocity*/
    fwf     = sf_output("out");/*wavfield snap*/
    frec    = sf_output("rec");/*record*/

    fG  = sf_input("G");
    fsx = sf_input("sx");
    fsz = sf_input("sz");
    
    if (SF_FLOAT != sf_gettype(fvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fsource)) sf_error("Need float input");

    /* parameters of souce */
    
    if (!sf_getbool("srcdecay", &srcdecay)) srcdecay=false;
    /* source decay */
    if (!sf_getint("srcrange", &srcreange)) srcrange=10;
    /* source decay range */
    if (!sf_getfloat("srctrunc", &srctrunc)) srctrunc=100;
    /* trunc source after srctrunc time (s) */
    
    if (!sf_getint("snapinter", &snapinter)) snapinter=1;
    /* snap interval */
    if (!sf_getbool("freesurface", &freesurface)) freesurface=false;
    
    /* Read/Write axes */
    at = sf_iaxa(fsource, 1); nt = sf_n(at); dt = sf_d(at);
    ax = sf_iaxs(fvel, 2);    nxb = sf_n(ax); dx = sf_d(ax); ox = sf_o(ax);
    az = sf_iaxs(fvel, 1);    nzb = sf_n(ax); dz = sf_d(az); oz = sf_o(az);

    if (!sf_histint(fGx, "n1", &nxz)) sf_error("No n1= in input");
    if (nxz != nxb*nzb) sf_error (" Need nxz = nxb*nzb");
    if (!sf_histint(fG, "n2", &LEN)) sf_error("No n2= in input");
    
    /*source loaction parameters*/
    slx = -1.0; spx = -1;
    slz = -1.0; spz = -1;
    gdep = -1.0; gp = 0;
    
    if (!sf_getfloat("slx", &slx)) ; 
    /*source location x */
    if (!sf_getint("spx", &spx));
    /*source location x (index)*/
    if((slx<0 && spx <0) || (slx>=0 && spx >=0 ))  sf_error("Need src location");
    if (slx >= 0 )    spx = (int)((slx-ox)/dx+0.5);
    
    if (!sf_getfloat("slz", &slz)) ;
    /* source location z */
    if (!sf_getint("spz", &spz)) ;
    /*source location z (index)*/
    if((slz<0 && spz <0) || (slz>=0 && spz >=0 ))  sf_error("Need src location");
    if (slz >= 0 )    spz = (int)((slz-ox)/dz+0.5);
    
    if (!sf_getfloat("gdep", &gdep)) ;
    /* recorder depth on grid*/
    if (!sf_getint("gp", &gp)) ;
    /* recorder depth on index*/
    if ( gdep>=oz) { gp = (int)((gdep-oz)/dz+0.5);}
    if (gp < 0.0) sf_error("gdep need to be >=oz");
    /*source and receiver location*/

    /* read FD schemes */
    sxtmp = sf_floatalloc(LEN);
    sztmp = sf_floatalloc(LEN);
    
    sx = sf_intalloc(LEN);
    sz = sf_intalloc(LEN);
    
    sf_floatread(sxtmp, LEN, fsx);
    sf_floatread(sztmp, LEN, fsz);
    
    G = sf_floatalloc3(nzb, nxb, LEN);
    
    sf_floatread(G[0][0], nzb*nxb*LEN, fG);

    mx = 0; mz = 0;
    for (ix=0; ix<LEN; ix++) {
	sx[ix] = (int)sxtmp[ix];
	mx = abs(sx[ix])>mx? abs(sx[ix]):mx;
    }
    for (iz=0; iz<LEN; iz++) {
	sz[iz] = (int)sztmp[iz];
	mz = abs(sz[iz])>mz? abs(sz[iz]):mz;
    }
    marg = mx>mz?mx:mz;

    nx = nxb -2*marg;
    nz = nzb -2*marg;

    free(sxtmp); free(sztmp);

    /* source and receiver location */
    sf_setn(ax, nx);
    sf_oaxa(frec, at, 1);
    sf_oaxa(frec, ax, 2);
    
    /* set axis for snap file */
    sf_setn(az, nz);
    sf_setn(at, (int)(nt-1)/snapinter+1);
    sf_setd(at, dt*snapinter);
    sf_oaxa(fwf, az, 1);
    sf_oaxa(fwf, ax, 2);
    sf_oaxa(fwf, at, 3);
    
    /* read model */
    vel = sf_floatalloc2(nzb, nxb);
    sf_floatread(vel[0], nxz, fvel);

    pn1 = sf_floatalloc2(nzb, nxb);
    pn0 = sf_floatalloc2(nzb, nxb);

    record = sf_floatalloc2(nt, nx);
    
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=0; ix<nxb; ix++) {
	for (iz=0; iz<nzb; iz++) {
	    pn0[ix][iz] = 0.0;
	}
    }

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=0; ix<nxb; ix++) {
	for (iz=0; iz<nzb; iz++) {
	    pn1[ix][iz] = 0.0;
	}
    }    




































