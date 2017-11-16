#ifndef gmx_nitrile_hbond_hpp
#define gmx_nitrile_hbond_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
#include <fstream>
#include <algorithm>


#ifdef __cplusplus
extern "C"
#endif
#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/statutil.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/trajana.h>

#ifndef MPI
#define MPI (4.*atan(1.))
#endif

#define RAD2DEG(X) (X*180./MPI)
#define DEG2RAD(X) (X*MPI/180.)

struct t_mol
{
    int resid;
    bool is_hb; 
    float nh;
    float cnh;
    float nho;
};

typedef std::vector<t_mol> t_water;
typedef std::vector<t_mol> t_prot ;

/*! \brief
 * Template analysis data structure.
 */
typedef struct
{
    float                               cr ; 
    float                               r  ; 
    float                               nho; 
    float                               cnh; 
} t_cutoffs ; 

typedef struct
{
    gmx_ana_selection_t                 *refsel;
    FILE                                *fp;
    FILE                                *fpp;
    FILE                                *fpa;
    FILE                                *fpg;
    FILE                                *fpnwg;                 //non-water geometry   

    int                                 framen;                 //number of frames 
    int                                 a1, a2;                 //indices of C and N, respectively, for nitrile
    int                                 nhb;                    //total number of hbonds 
    int                                 nwater;                 //number of hbonds to water
    int                                 nprot ;                 //number of hbonds to protein 

    std::vector<int>                    p_water;                //number of persistent waters
    std::vector<int>                    p_hbonds;               //number of persistent hbonds from water 
    std::vector<int>                    p_prot_hb;              //number of persistent hbonds from protein

    std::vector<t_water>                water;                  //vector of mol objects for each nearby water 
    std::vector<t_water>                water_hb;               //vector of mol objects for each hbonding water
    std::vector<t_prot>                 prot ;                  //vector of mol objects for each hbonding protein 

    std::vector<int>                    frame_sigma_nhb;        //number of Cho sigma water hbonds for each frame
    std::vector<int>                    frame_pi_nhb;           //number of Cho pi water hbonds for each frame
    std::vector<int>                    frame_lQ_nhb;           //number of LQ water hbonds for each frame  
    std::vector<int>                    frame_prot_nhb;         //number of protein water hbonds for each frame 
    std::vector<int>                    frame_totnhb;           //number of all hbonds for each frame 

    t_cutoffs                           cutoffs ; 
    bool                                bVerbose;
    bool                                doPersistent, doLog, doWatGeo,doProtGeo;
} t_analysisdata;

bool parse_water(const float *a1, const float *a2, const float *ow, const float *h1, const float* h2, t_cutoffs cutoffs, t_mol &mol);

bool parse_Hs(const float *a1, const float *a2, const float *ow, const float *h1, t_cutoffs cutoffsr, t_mol &mol);

bool parse_Hs(const float *a1, const float *a2, const float *ow, const float *h1, const float *h2, t_cutoffs cutoffs, t_mol &mol);

bool parse_Hs(const float *a1, const float *a2, const float *ow, const float *h1, const float *h2, const float *h3, t_cutoffs cutoffs, t_mol &mol);

int count_persistant(std::vector<std::vector<t_mol> > mol, int framen, std::vector<int> &persistant) ; 

void vsub( const float *a1, const float *a2, const int &size, float *d );

float vlen_( const float *v );

float vlen(const float *a1, const float *a2 );

float vangle(const float *a1, const float *a2, const float *a3);

float ndot( const float *a1, const float *a2 );

#endif
