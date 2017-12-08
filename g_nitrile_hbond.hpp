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
    //int resid;     // internal Gromacs residue index 
    int resNumber; // residue Number that is printed on topolgy, VMD, .gro files, etc. 
    bool is_hb;    // is it hydrogen bonding
    float nh;      // Donor - Hydrogen distance   (N:--H) 
    float cnh;     // X - Acceptor - Hydrogen i   (C=N:--H) 
    float nho;     // Acceptor - Hydrogen - Donor (N:--H-O) 
};

typedef std::vector<t_mol> t_water;
typedef std::vector<t_mol> t_prot ;

/*! \brief
 * Template analysis data structure.
 */

typedef struct
{
    float                               cr ;                    //Cutoff for HBond donor to be considered close to hbond. 
                                                                //      If within this distance it will be counted a 
                                                                //      Cho-sigma or Cho-pi. No angle requirements. 
    float                               r  ;                    //Cutoff for LeQuestial Hydrogen bond. Also depends on angle requirements below. 
    float                               cnh;                    //  X      - Acceptor - Hydrogen angle (C=N:--H) cutoff for LQ HBond. 
    float                               nho;                    //Acceptor - Hydrogen - Donor    angle (N:--H-O) cutoff for LQ HBond. 
} t_cutoffs ; 

typedef struct
{
    gmx_ana_selection_t                 *refsel;                //GMX struct that contains the structural information of a Group in a particular frame
    FILE                                *fp;                    //File (write) pointer to hbonding per frame information
    FILE                                *fpp;                   //File (write) pointer to hbonding persistancy data
    FILE                                *fpa;                   //File (write) pointer to accumulated hbonding
    FILE                                *fpg;                   //File (write) pointer to water hbonding geometry. 
    FILE                                *fpnwg;                 //File (write) pointer to non-water hbonding geometry   

    int                                 framen;                 //counter to count the number of frames analyzed. 
    int                                 a1, a2;                 //indices of C and N, respectively, for nitrile
    int                                 nhb;                    //counter of total number of hbonds 
    int                                 nwater;                 //counter of number of hbonds to water
    int                                 nprot ;                 //counter of number of hbonds to protein 

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

    t_cutoffs                           cutoffs ;               //Struct to hold all cutoff values 
    bool                                bVerbose;               //Not impolemented yet. 
    bool                                doPersistent, doLog, doWatGeo,doProtGeo;
} t_analysisdata;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Considers a molecule with one possible hydrogen, calculates if the molecule is hydrogen                                     // 
//            bonding or not, and sets the geometry properties of the molecule                                                            //
//Arguments:  Accepts points to arrays of the coordinates of the C of the nitrile, the N of the nitrile, the Donor (O in water)           //
//            and the Hydrogen respectively, then the struct of cutoff information, and a reference to the struct of molecule information.// 
//Return:     Returns true if the molecule is "close" enough to be considered Cho bonding, and false if outside the Cho radius.           // 
//            The mol struct is passed by reference and has it's is_hb, cnh, nw, nho propterties calculated and modified.                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool parse_Hs(const float *a1, const float *a2, const float *ow, 
              const float *h1, 
              t_cutoffs cutoffsr, t_mol &mol);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Considers a molecule with two possible hydrogens, and using the nearest hydrogen, calculates if the molecule is hydrogen    // 
//            bonding or not, and sets the geometry properties of the molecule                                                            //
//Arguments:  Accepts points to arrays of the coordinates of the C of the nitrile, the N of the nitrile, the Donor (O in water)           //
//            and the two Hydrogens respectively, then the struct of cutoff information,                                                  //
//            and a reference to the struct of molecule information.                                                                      // 
//Return:     Returns true if the molecule is "close" enough to be considered Cho bonding, and false if outside the Cho radius.           // 
//            The mol struct is passed by reference and has it's is_hb, cnh, nw, nho propterties calculated and modified.                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool parse_Hs(const float *a1, const float *a2, const float *ow, 
              const float *h1, const float *h2, 
              t_cutoffs cutoffs, t_mol &mol);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Considers a molecule with three possible hydrogens, and using the nearest hydrogen, calculates if the molecule is hydrogen  // 
//            bonding or not, and sets the geometry properties of the molecule                                                            //
//Arguments:  Accepts points to arrays of the coordinates of the C of the nitrile, the N of the nitrile, the Donor (O in water)           //
//            and the two Hydrogens respectively, then the struct of cutoff information,                                                  //
//            and a reference to the struct of molecule information.                                                                      // 
//Return:     Returns true if the molecule is "close" enough to be considered Cho bonding, and false if outside the Cho radius.           // 
//            The mol struct is passed by reference and has it's is_hb, cnh, nw, nho propterties calculated and modified.                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool parse_Hs(const float *a1, const float *a2, const float *ow, 
              const float *h1, const float *h2, const float *h3, 
              t_cutoffs cutoffs, t_mol &mol);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Considers an array of molecules, and counts number of times a molecule appears in consecutive frames (persistantly).        //
//            Count the frequency of a molecule appearing persistantly as a function of time.                                             //
//Arguents:   Accepts an array of molecule structs, and integer of how many frames are available, and an integer vector by reference.     //
//Returns :   Returns 0 if successful. Integer array by reference is modified to contain the number of times a molecule was persistent    //
//            for that amount of time, at each index of the array                                                                         //  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_persistant(std::vector<std::vector<t_mol> > mol, int framen, std::vector<int> &persistant) ; 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Calculates the element by element difference of two vectors of arbitrary length.                                            //
//Arguments:  Pointer to two vectors, the size of both vectors (must be equal size), pointer to vector where the difference will stored.  //
//Returns:    Void. Pointer to difference vector modified so the difference can be returned.                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void vsub( const float *a1, const float *a2, const int &size, float *d );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Calculates the magnitude of a vector of length 3.                                                                           //
//Arguments:  Pointer to a vector. Must be of length 3! (x,y,z vector)                                                                    //
//Returns:    Magnitude of vector as a float.                                                                                             //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float vlen_( const float *v );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Calculates the distance between two three coordinate vectors.                                                               //
//Arguments:  Pointer to two vectors. Must be of length 3! (x,y,z vector)                                                                 //
//Returns:    Distance between two vectors as a float.                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float vlen(const float *a1, const float *a2 );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Calculates the angle between three points in space.                                                                         //
//Arguments:  Three pointers to coordinates of each point in space (x,y,z vector)                                                         //
//Returns:    The angle between the three points as a float.                                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float vangle(const float *a1, const float *a2, const float *a3);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Purpose:    Calculates the dot product between two vectors of length three.                                                             //
//Arguments:  Two pointers to three-coordinate vectors.                                                                                   // 
//Returns:    The dot product of the two vectors as a float.                                                                              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float ndot( const float *a1, const float *a2 );

#endif
