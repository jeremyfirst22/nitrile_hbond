/*
 * This file is part of the GROMACS molecular simulation package.
 ource /usr/local/gromacs/bin/GMXRC
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009,2010,2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "g_nitrile_hbond.hpp"

/*! \brief
 * Function that does the analysis for a single frame.
 *
 * It is called once for each frame.
 */
static int analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata      *d = (t_analysisdata *)data;
    t_water             water;
    t_prot              prot ; 
    t_water             water_hb ; 
    int sigma_nhb = 0, pi_nhb = 0, lQ_nhb = 0, prot_nhb = 0, totnhb = 0; 
    float a1[3] = { fr->x[d->a1][0], fr->x[d->a1][1], fr->x[d->a1][2] };
    float a2[3] = { fr->x[d->a2][0], fr->x[d->a2][1], fr->x[d->a2][2] };
    int atom_ndx ; 
    int resid ; 
    int resNumber ; 
    char* resname ; 
    char* atomname ; 
    t_mol mol;                                        //create mol object 
    // Loop through everything   
    for (int g = 0; g < nr; ++g) {                    // for each group
        if (g > 1) {
            fprintf(stderr,"\nwarning: more than 1 group was specified.  cowardly ignoring all additional groups\n");
            break;
        }
        for (int i = 0; i < sel[g]->p.nr; ++i) {      // for each atom in the group
            mol.is_hb = false ; 
            //mol.resid = resid ; 

            atom_ndx = sel[g]->g->index[i];           // how we get to the atom index
            resid = top->atoms.atom[atom_ndx].resind; // internal residue index 
            resNumber = top->atoms.resinfo[resid].nr ;
            resname = *top->atoms.resinfo[resid].name;
            atomname = *top->atoms.atomname[atom_ndx];
            char elem = atomname[0] ; 

            mol.resNumber = resNumber ; 

            //Water molecules must be either of SOL of HOH 
            if (strncmp(resname, "SOL", 4) == 0 || strncmp(resname, "HOH", 4) == 0) { 
                if (strncmp(atomname, "OW", 3) == 0) { 
                    //fprintf(stdout, "\nresid: %i resname: %s atomname: %s atom_ndx: %i\n", resid, resname, atomname, atom_ndx) ; 
                    //Get coords of water atoms
                    float OW[3] = { fr->x[atom_ndx][0], fr->x[atom_ndx][1], fr->x[atom_ndx][2] };
                    float H1[3] = { fr->x[atom_ndx+1][0], fr->x[atom_ndx+1][1], fr->x[atom_ndx+1][2] };
                    float H2[3] = { fr->x[atom_ndx+2][0], fr->x[atom_ndx+2][1], fr->x[atom_ndx+2][2] };

                    //parse_water returns true if water is close enough to possible h-bond
                    if (parse_Hs(a1, a2, OW, H1, H2, d->cutoffs, mol)) { 
                        if (mol.is_hb) { 
                            totnhb++;                 //counter for all hydrogen bonds
                            lQ_nhb++;                 //specifically water hydrogen bonds
                            water_hb.push_back(mol) ; //add to water hbond array
                            //fprintf(stdout, "\n\tFrame: %i  Water %i is hbonding\n",d->framen, atom_ndx) ; 
                            //fprintf(stdout, "\n\tFrame: %i  Water %i is hbonding\n",d->framen, resid) ; 
                            //fprintf(stdout, "\n\tFrame: %i  Water %i is hbonding\n",d->framen, resNumber) ; 
                        }
                        if (mol.cnh >= 120.) {        //Cho hbonds - only a distance requirement 
                            sigma_nhb++;              //Decide if these are sigma or pi 
                        }                             //Note: These do NOT count as actual hbonds, 
                        else {                        //      this is just to compare against Cho paper
                            pi_nhb++;
                        }
                        water.push_back(mol);
                        d->nwater++;
                    }
                    //fprintf(stdout, "\tAdded water. length of water = %i\n",water.size())  ; 
                }
            }
            else { 
                //Any N, O, S, or CG on Met are the allowed donors as of now. 
                if (elem == 'N' || elem == 'O' || elem == 'S' || strncmp(atomname, "CA",3) == 0) { 
    //                fprintf(stdout,"%s found in group! \t%s\t%i\n",atomname, resname, resNumber) ; 
                    //fprintf(stdout, "\t\t atom_ndx = %i\n",atom_ndx) ; 

                    //Coords of Hbond donor  
                    float donor[3] = { fr->x[atom_ndx][0], fr->x[atom_ndx][1], fr->x[atom_ndx][2] } ;
                    int maxHs = 3  ;            //Maximum Hs attached to the donors above is 3. 
                    int countHs = 0 ; 
                    for (countHs = 0 ; countHs < maxHs ; countHs++) {
                        char *nextAtom = *top->atoms.atomname[atom_ndx+countHs+1] ; 
                        char c = nextAtom[0] ; 
    //                    fprintf(stdout, "\t\tResidue: %s\tAtomname: %s\tElement: %c\n",resname, nextAtom,c) ; 
                        if ( c != 'H' ) {
                            //fprintf(stdout, "\t\t\tNot hydrogen, breaking loop\n") ; 
                            break ; 
                        }
                    }
                    //fprintf(stdout, "%i\n",countHs) ; 

                    //These are the xyz coords of the next three atoms, but they may or may not be Hs
                    float H1[3] = { fr->x[atom_ndx+1][0], fr->x[atom_ndx+1][1], fr->x[atom_ndx+1][2] } ;
                    float H2[3] = { fr->x[atom_ndx+2][0], fr->x[atom_ndx+2][1], fr->x[atom_ndx+2][2] } ;
                    float H3[3] = { fr->x[atom_ndx+3][0], fr->x[atom_ndx+3][1], fr->x[atom_ndx+3][2] } ; 
                    switch (countHs){
                        case 0 : 
                            //Not an H donor since no Hs attached. 
    //                        fprintf(stdout,"\t\t\tCase 0\n") ; 
                            break ; 
                        case 1 : 
                            parse_Hs(a1, a2, donor, H1, d->cutoffs, mol) ;  
    //                        fprintf(stdout,"\t\t\tCase 1\n") ; 
                            break ; 
                        case 2: 
                            parse_Hs(a1, a2, donor, H1, H2, d->cutoffs, mol) ;  
    //                        fprintf(stdout,"\t\t\tCase 2\n") ; 
                            break ; 
                        case 3: 
                            parse_Hs(a1, a2, donor, H1, H2, H3, d->cutoffs, mol) ;  
    //                        fprintf(stdout,"\t\t\tCase 3\n") ; 
                            break ; 
                    }

                    //fprintf(stdout, "\t\t\tNH= %6.4f %6.4f %6.4f\n",a2[0],a2[1],a2[2]) ; 
                    //fprintf(stdout, "\t\t\tN = %6.4f %6.4f %6.4f\n",N[0],N[1],N[2]) ; 
                    //fprintf(stdout, "\t\t\tH = %6.4f %6.4f %6.4f\n",H[0],H[1],H[2]) ; 
                    //fprintf(stdout, "\tCalling parse_nonwater\n") ; 
                    if (mol.is_hb) { 
                        //if (d->doGeo) {
                        //    fprintf(d->fpg, "%10i %10i %12.4f %12.4f %12.4f\n", d->framen, atom_ndx, 10*mol.nh, mol.cnh, mol.nho);
                        //}
                        totnhb++;                 //counter for all hydrogen bonds
                        prot_nhb++ ;              //specifically protein hydrogen bonds
                        prot.push_back(mol) ;     //add to protein hbond array
                        //fprintf(stdout, "\n\tFrame: %i  Residue %i is hbonding. Donor is %s\n",d->framen, resNumber,atomname) ; 
                    }
                }
            }
        }
    }
    d->water.push_back(water);
    d->water_hb.push_back(water_hb);
    d->prot.push_back(prot);
    //fprintf(stdout, "\nFrame complete. prot.size() = %i\n",prot.size() ) ; 
    d->frame_lQ_nhb.push_back(lQ_nhb);
    d->frame_sigma_nhb.push_back(sigma_nhb);
    d->frame_pi_nhb.push_back(pi_nhb);
    d->frame_prot_nhb.push_back(prot_nhb) ;
    d->frame_totnhb.push_back(totnhb) ;
    //fprintf(stdout, "\n\t\tFrame = %i Total hbonds = %i Prot = %i \n",d->framen,totnhb,prot_nhb) ; 
    //fprintf(stdout, "\t\t\tLQ Water = %i Sigma = %i Pi = %i\n\n", lQ_nhb,sigma_nhb, pi_nhb) ; 
    /* increment the frame number */
    d->framen++;
    /* We need to return 0 to tell that everything went OK */
    return 0;
}

int analyze_information(void *data) {
    t_analysisdata      *d = (t_analysisdata *)data;
    int maxhb = 6;
    maxhb++ ;       //It does not make sense to count from zero. 

    //We will use these vectors to count through the number of frames with each number of Hbonds
    std::vector<int> frames_with_nearbywater(maxhb,0);
    std::vector<int> frames_with_totnhb (maxhb,0);
    std::vector<int> frames_with_lQ (maxhb,0);
    std::vector<int> frames_with_sigma (maxhb,0);
    std::vector<int> frames_with_pi (maxhb,0);
    std::vector<int> frames_with_prot (maxhb,0);

    //We will use these vectors to accumulate the persistant hbonds
    std::vector<int> persistant_water (d->framen,0);
    std::vector<int> persistant_hb (d->framen,0);
    std::vector<int> persistant_prothb (d->framen,0);

    //We will use these vectors to count lifetimes of hbonds
    std::vector<int> total_water (d->framen,0);
    std::vector<int> total_hb (d->framen,0);
    std::vector<int> total_prothb (d->framen,0);

    /*fprintf(stdout, "\nPersist waters: %i  ",persistant_water.size() ) ; 
    for (int i = 0 ; i < persistant_water.size() ; i++){
        fprintf(stdout, "%i  ",persistant_water[i]) ; 
    }*/
    /*fprintf(stdout, "\nPersist hb:     %i  ",persistant_hb.size()) ; 
    for (int i = 0 ; i < persistant_hb.size() ; i++){
        //fprintf(stdout, "%i  ",persistant_hb[i]) ; 
    }*/
    /*fprintf(stdout, "\nPersist prot:   %i  ",persistant_prothb.size()) ;
    for (int i = 0 ; i < persistant_prothb.size() ; i++) {
        //fprintf(stdout, "%i  ",persistant_prothb[i]) ; 
    }*/

    //Print number hbond information for each frame to a file
    //    Default = frame_hb.xvg 
    for (int i=0; i<d->framen; i++) {
        fprintf(d->fp,"%10i %10i %10i %10i %10i %10i %10i\n",i,d->water[i].size(), d->frame_lQ_nhb[i], d->frame_sigma_nhb[i], d->frame_pi_nhb[i], d->frame_prot_nhb[i], d->frame_totnhb[i]); 
    }  

    //Analyze number of frames with different number of hbonds
    //   Default = hb_count.xvg  
    if (d->doLog) {
        //Count frames with different number of hbonds. 
        for (int i = 0; i<d->framen; i++) {
            for (int j=0; j<maxhb; j++) {
                if (d->water[i].size() == j) {
                    frames_with_nearbywater[j] ++; 
                }
                if (d->frame_totnhb[i] == j) {
                    frames_with_totnhb[j] ++; 
                }
                if (d->frame_lQ_nhb[i] == j) {
                    frames_with_lQ[j] ++;
                }
                if (d->frame_sigma_nhb[i] == j) {
                    frames_with_sigma[j] ++;
                }
                if (d->frame_pi_nhb[i] == j) {
                    frames_with_pi[j] ++;
                }
                if (d->frame_prot_nhb[i] == j) {
                    frames_with_prot[j] ++;
                }
            }
        }
        //Print to file. 
        for (int i=0;i<maxhb;i++) {
            fprintf(d->fpa,"%10i %10i %10i %10i %10i %10i %10i\n",i,frames_with_nearbywater[i],frames_with_lQ[i],frames_with_sigma[i],frames_with_pi[i],frames_with_prot[i],frames_with_totnhb[i]);
        }
    }

    //Analyze Geometry of each hydrogen-bonding water molecule 
    //   Default = wat_geometry.xvg 
    if (d->doWatGeo) {
        //For each frame
        for (int i = 0 ; i < d->framen ; i++) { 
            //For each molecule in each frame
            for (int j = 0 ; j < d->water_hb[i].size() ; j++) {
                //Print to file 
                fprintf(d->fpg, "%10i %10i %12.4f %12.4f %12.4f\n",i,d->water_hb[i][j].resNumber,10*d->water_hb[i][j].nh,d->water_hb[i][j].cnh,d->water_hb[i][j].nho) ; 
            }
        }
    }

    //Analyze Geometry of each hydrogen-bonding protein residue  
    //   Default = prot_geometry.xvg 
    if (d->doProtGeo) {
        //For each frame
        for (int i = 0 ; i < d->framen ; i++) { 
            //For each molecule in each frame
            for (int j = 0 ; j < d->prot[i].size() ; j++) {
                //Print to file 
                fprintf(d->fpnwg, "%10i %10i %12.4f %12.4f %12.4f\n",i,d->prot[i][j].resNumber,10*d->prot[i][j].nh,d->prot[i][j].cnh,d->prot[i][j].nho) ; 
            }
        }
    }

    if (d->doPersistent) {
        //fprintf(stdout, "\n\n\t\t***Now getting persistent water information.***\n") ; 
        //fprintf(stdout, "\n\t\t\t\td->framen = %i\n\n",d->framen) ; 

        count_persistant(d->water, d->framen, persistant_water,d->forgivenessLevel) ; 
        count_persistant(d->water_hb,d->framen,persistant_hb,d->forgivenessLevel) ; 
        count_persistant(d->prot,d->framen,persistant_prothb,d->forgivenessLevel) ; 
          
        //Accumulate how lifetimes
        for (int i=0; i<d->framen; i++) {
            for (int j=i; j<d->framen; j++) {
                total_water[i] += persistant_water[j];
                total_hb[i] += persistant_hb[j];
                total_prothb[i] += persistant_prothb[j] ; 
            }
        // Write to xvg file   
        for (int i=0; i<d->framen;i++){
            fprintf(d->fpp,"%10i %10i %10i %10i %10i %10i %10i\n",i,persistant_water[i],persistant_hb[i],persistant_prothb[i],total_water[i], total_hb[i],total_prothb[i]);
            }
        }
    }

    /*fprintf(stdout, "\nPersist waters: %i  ",persistant_water.size() ) ; 
    for (int i = 0 ; i < persistant_water.size() ; i++){
        fprintf(stdout, "%i  ",persistant_water[i]) ; 
    }*/
    /*fprintf(stdout, "\nPersist hb:     %i  ",persistant_hb.size()) ; 
    for (int i = 0 ; i < persistant_hb.size() ; i++){
        //fprintf(stdout, "%i  ",persistant_hb[i]) ; 
    }*/
    /*fprintf(stdout, "\nPersist prot:   %i  ",persistant_prothb.size()) ;
    for (int i = 0 ; i < persistant_prothb.size() ; i++) {
        //fprintf(stdout, "%i  ",persistant_prothb[i]) ; 
    }*/
    return 0 ; 
}

int count_persistant(std::vector<std::vector<t_mol> > mol, int framen, std::vector<int>& time_persistant, int forgivenessLevel){
    //fprintf(stdout, "\n***Array size = %i***\n\n",mol.size() ) ; 
    for (int i = 0 ; i < framen ; i++) { 
        //fprintf(stdout, "Time = %i\n",i) ; 
        //Loop through hydrogen bonding waters in each frame 
        for (int j = 0 ; j<mol[i].size() ; j++){
            //Find out if mol was in previous frame
            bool previous = false ; 
            //Loop through all mols in previous frame to see if this mol was in last frame
            if (i>0) { //If first frame, must be new mol
                int f = 0 ; 
                do {
                    //fprintf(stdout, "\tSearching %i frames back (k = %i)\n", f+1, i-f-1) ; 
                    //fprintf(stdout, "\t\tChecking %i molecules: ",mol[i-f-1].size() ) ; 
                    for (int k = 0 ; k < mol[i-f-1].size() ; k++) {
                        //fprintf(stdout, "%i... ",k+1) ; 
                        if(mol[i-f-1][k].resNumber == mol[i][j].resNumber ) {
                            previous = true ; 
                            //fprintf(stdout, "\t\tFound %i frames back\n",f+1) ; 
                            break ; 
                        }
                    }
                    //fprintf(stdout, "\n") ; 
                    f++ ; 
                } while (f <= forgivenessLevel && !previous && i-f-1 > 0) ; 
            }
            //If a new residue, then find out how long it lasts
            if (! previous) { 
                //fprintf(stdout, "\tNew residue: %i\n", mol[i][j].resNumber) ; 
                int persFrames = 0 ; //counter for number of persistent frames

                //fprintf(stdout, "\t\tChecking:\n") ; 
                bool persistant = true  ; 
                int k = 0 ; 
                while (persistant && i+k+1 < framen ) { 
                    //fprintf(stdout, "\t\t    t=%i ",i+1+k) ; 
                    persistant = false ; 
                    for (int l = 0 ; l < mol[i+k+1].size() ; l++) {
                        //fprintf(stdout, " k= %i res.= %i ", k, mol[k][l].resNumber) ; 
                        if(mol[i+k+1][l].resNumber == mol[i][j].resNumber ) {
                            //fprintf(stdout, " Found\n")  ; 
                            persistant = true ; 
                            persFrames++ ; 
                            break ; 
                        }
                    }
                    if ( ! persistant ) {
                        //fprintf(stdout, "not found\n") ; 
                        for (int f = 1 ; f <= forgivenessLevel && i+k+f+1 < framen; f++) {
                            //fprintf(stdout, "\t\t\tChecking f level: %i ",f) ;  
                            for (int l = 0 ; l < mol[i+k+f+1].size() ; l++) {
                                if(mol[i+k+f+1][l].resNumber == mol[i][j].resNumber ) {
                                    //fprintf(stdout, "found\n") ; 
                                    persistant = true ; 
                                    persFrames++ ; 
                                    break ; 
                                }
                            }
                            if (persistant) break ; 
                            //else fprintf(stdout, "not found\n") ; 
                        } 
                        if(!persistant ) {
                            //fprintf(stdout, "\t\t\tOut of tries. Not persistant.\n") ; 
                        }
                    } 
                    k++ ; 
                }   
                //fprintf(stdout, "\n") ; 
                //fprintf(stdout, "\t\t\tPersistant for %i frames\n",persFrames) ; 
                time_persistant[persFrames]++ ; 
            }
        } //end mols
    }
    return 0 ; 
}

bool parse_Hs(const float *a1, const float *a2, const float *donor, 
    const float *h1, 
    t_cutoffs cutoffs, t_mol &mol)
{
    float v1 = vlen(a2,h1);
    if (v1 <= cutoffs.cr ) { 
        //Only one hydrogen - donor length to consider
        mol.nh = v1;
        mol.cnh = vangle(a1,a2,h1);
        mol.nho = vangle(a2,h1,donor);
        if (mol.nh < cutoffs.r && mol.cnh > cutoffs.cnh && mol.nho > cutoffs.nho ) { 
            mol.is_hb++ ;           //true 
        } 
        return true;
    }
    return false;
}

bool parse_Hs(const float *a1, const float *a2, const float *donor, 
    const float *h1, const float *h2, 
    t_cutoffs cutoffs, t_mol &mol)
{
    float v1 = vlen(a2,h1);
    float v2 = vlen(a2,h2);

    if (v1 <= cutoffs.r || v2 <= cutoffs.r) 
    {
        if (v1 < v2) {             //v1 is the smallest 
            mol.nh = v1;
            mol.cnh = vangle(a1,a2,h1);
            mol.nho = vangle(a2,h1,donor);
        }
        else {                     //v2 is the smallest
            mol.nh = v2;
            mol.cnh = vangle(a1,a2,h2);
            mol.nho = vangle(a2,h2,donor);
        }
        if (mol.nh < cutoffs.r && mol.cnh > cutoffs.cnh && mol.nho > cutoffs.nho ) { 
            mol.is_hb++ ;          //true 
        } 
        return true;
    }
    return false;
}

bool parse_Hs(const float *a1, const float *a2, const float *donor, 
    const float *h1, const float *h2, const float *h3, 
    t_cutoffs cutoffs, t_mol &mol)
{
    float v1 = vlen(a2,h1);
    float v2 = vlen(a2,h2);
    float v3 = vlen(a2,h3);

    if (v1 <= cutoffs.cr || v2 <= cutoffs.cr || v3 <= cutoffs.cr) { 
        if (v1 < v2 && v1 < v3) {  //v1 is the smallest 
            mol.nh = v1;
            mol.cnh = vangle(a1,a2,h1);
            mol.nho = vangle(a2,h1,donor);
        }
        else if (v2 < v3) {        //v2 is the smallest
            mol.nh = v2;
            mol.cnh = vangle(a1,a2,h2);
            mol.nho = vangle(a2,h2,donor);
        }
        else {                     //v3 is the smallest
            mol.nh = v2;
            mol.cnh = vangle(a1,a2,h2);
            mol.nho = vangle(a2,h2,donor);
        }
        if (mol.nh < cutoffs.r && mol.cnh > cutoffs.cnh && mol.nho > cutoffs.nho) {
            mol.is_hb++ ;          //true 
        } 
        return true;
    }
    return false;
}

void vsub( const float *a1, const float *a2, const int &size, float *d )
{
    for (int i=0; i<size; i++) {
        d[i] = a2[i] - a1[i];
    }
}

float vlen_( const float *v )
{
    return pow( v[0]*v[0] + v[1]*v[1] + v[2]*v[2], 0.5);
}

float vlen(const float *a1, const float *a2 )
{
    float vec[3];
    vsub(a1,a2,3,vec);
    return vlen_(vec);
}

float ndot( const float *a1, const float *a2 )
{
    float rv1 = 1./vlen_(a1);
    float rv2 = 1./vlen_(a2);
    return a1[0]*rv1*a2[0]*rv2 + a1[1]*rv1*a2[1]*rv2 + a1[2]*rv1*a2[2]*rv2;
}

float vangle(const float *a1, const float *a2, const float *a3)
{
    float v21[3], v23[3];
    vsub(a2,a1,3,v21);
    vsub(a2,a3,3,v23);
    float d = ndot(v21,v23);
    if (d > 1) { d = 1; };
    return RAD2DEG(acos(d));
}

/*! \brief
 * Function that implements the analysis tool.
 *
 * Following the style of Gromacs analysis tools, this function is called
 * \p gmx_something.
 */

int gmx_nitrile_hbond(int argc, char *argv[])
{
    const char *desc[] = {
        "\tComputes the number of hydrogen bonds to a hydrogen bond acceptor",
        "built specifically for nitrile labels in proteins",
        "Use the -a1 flag for the hydrogen bond acceptor (C)",
        "and the -a2 flag for the electronegative atom (N).",
        "\n\tUse the -select flag to specify which atoms you",
        "want to parse through.  This can have a large effect",
        "on performance, and as such should be used.",
        "Using -select",
        "'resname SOL and same residue as within 0.5 of resname",
        "CNC and name NE' is recommended because it would loop over",
        "only waters near the cyanocysteine nitrogen.",
        "\n"
    };
    /* Command-line arguments */
    /*
     Hydrogen-bond acceptor properties of nitriles: a combined crystallographic and ab initio theoretical investigation
     Jean-Yves Le Questel, Michel Berthelot and Christian Laurence
     r average = .205(.020)
     theta1 = CNH, average = 145(23)
     theta2 = NHO, average = 156(18)
     95% of data within +/- 2 * STD
     */
    t_analysisdata      d;
    d.a1        = -1;
    d.a2        = -1;
    d.forgivenessLevel = 0 ; 
    d.bVerbose  = false;
    d.cutoffs.cr        = 0.30;
    d.cutoffs.r         = 0.205 + (2*0.02);
    d.cutoffs.cnh       = 145 - (2*23);
    d.cutoffs.nho       = 156 - (2*18);
    d.nhb       = 0;
    
    t_pargs         pa[] = {
        { "-a1", TRUE, etINT,
            {&d.a1}, "Starting atom for bond vector--ie: CD in CNC"},
        { "-a2", TRUE, etINT,
            {&d.a2}, "Ending atom for bond vector--ie: NE in CNC"},
        { "-forgiveness", FALSE, etINT, 
            {&d.forgivenessLevel}, "Forgiveness level for persistancy of hydrogen bonds (number of frames)"},
        { "-NHO_cutoff", FALSE, etREAL,
            {&d.cutoffs.nho}, "Minimum N-H-O bond angle to count as a hydorgen bond (degrees)"},
        { "-CNH_cutoff", FALSE, etREAL,
            {&d.cutoffs.cnh}, "Minimum N-H-O bond angle to count as a hydorgen bond (degrees)"},
        { "-angle_NH_cutoff", FALSE, etREAL,
            {&d.cutoffs.r}, "Minimum N-H bond distance to count as a hydrogen bond include angle criteria (nm)"},
        { "-distance_NH_cutoff", FALSE, etREAL,
            {&d.cutoffs.cr}, "Minimum N-H bond distance to count as a hydrogen bond with only distance criteria (nm)"},
        { "-v", FALSE, etBOOL,
            {&d.bVerbose}, "Be slightly more verbose"}
    };
    
    t_filenm        fnm[] = {
         { efXVG, "-o",  "frame_hb", ffWRITE }
        ,{ efXVG, "-op", "persistent", ffWRITE }
        ,{ efXVG, "-oa", "hb_count", ffWRITE }
        ,{ efXVG, "-or", "geometry", ffWRITE }
        ,{ efXVG, "-onwr", "nw_geometry", ffWRITE }
    };
    d.doLog = false, d.doPersistent = false, d.doWatGeo = false, d.doProtGeo = false ;
    
#define NFILE asize(fnm)
    
    gmx_ana_traj_t      *trj;
    t_topology          *top;
    output_env_t        oenv;
    int                 ngrps;
    gmx_ana_selection_t **sel;
    int                 g;
    int                 rc;
    
    CopyRight(stderr, argv[0]);
    
    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nanagrps(trj, -1);
    parse_trjana_args(trj, &argc, argv, 0,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);
    gmx_ana_get_topology(trj, FALSE, &top, NULL);

    gmx_ana_get_nanagrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);
    gmx_ana_init_coverfrac(trj, CFRAC_SOLIDANGLE);
    
    /* open xvg files */
    d.fp = NULL;
    d.fp = xvgropen(opt2fn("-o", NFILE, fnm), "Trajectory Hydrogen Bonds","Step", "Number of Water Molecules", oenv);
    const char * legend[] = {   "Frame Number",
                                "Nearby Water Molecules",
                                "Le Questel HBonded Water Molecules",
                                "Cho Sigma HBonded Water Molecules",
                                "Cho Pi HBonded Water Molecules", 
                                "Protein hydrogen bonds", 
                                "Total hydrogen bonds (LQ and protein)"
                            };
    xvgr_legend(d.fp,asize(legend),legend,oenv);

    d.fpa = NULL;
    if (opt2bSet("-oa", NFILE, fnm))
    {
        d.doLog = true;
        const char *alegend[] = {   "Number of each type per frame",
                                    "Number of Frames (Nearby Water Molecules)",
                                    "Number of Frames (Le Questel HBonded Water Molecules)",
                                    "Number of Frames (Cho Sigma HBonded Water Molecules)",
                                    "Number of Frames (Cho Pi HBonded Water Molecules)", 
                                    "Number of Frames (Protein hydrogen bonds)", 
                                    "Number of Frames (Total hydrogen bonds (LQ and protein))"
                                };
        d.fpa = xvgropen(opt2fn("-oa", NFILE, fnm), "Number of Frames with different numbers of Hydrogen Bonds", "Number of frames", "Number of Frames", oenv);
        xvgr_legend(d.fpa,asize(alegend),alegend,oenv);
    }

    d.fpg = NULL;
    if (opt2bSet("-or", NFILE, fnm))
    {
        d.doWatGeo = true;
        const char *glegend[] = {   "Frame",
                                    "Water Index",
                                    "NH Distance (Angstrom)",
                                    "CNH Angle (Degrees)",
                                    "NHO Angle (Degrees)"
                                };
        d.fpg = xvgropen(opt2fn("-or", NFILE, fnm), "Hydrogen Bonding Geometry", "Step", "Geometry", oenv);
        xvgr_legend(d.fpg,asize(glegend), glegend, oenv);
    }

    d.fpnwg = NULL;
    if (opt2bSet("-onwr", NFILE, fnm))
    {
        d.doProtGeo = true;
        const char *nwglegend[] = { "Frame",
                                    "Residue Index",
                                    "NH Distance (Angstrom)",
                                    "X - Acceptor - H Angle (CNH Angle) (Degrees)",
                                    "Acceptor - H - Donor Angle (Degrees) " 
                                };
        d.fpnwg = xvgropen(opt2fn("-onwr", NFILE, fnm), "Hydrogen Bonding Geometry", "Step", "Geometry", oenv);
        xvgr_legend(d.fpnwg,asize(nwglegend), nwglegend, oenv);
    }
    
    d.fpp = NULL;
    if (opt2bSet("-op", NFILE, fnm))
    {
        d.doPersistent = true;
        const char *flegend[] = {   "Persistent Water Molecules",
                                    "Persistent Hydrogen Bonds",
                                    "Persistent Hydrogen Bonds From Protein",
                                    "Total Persistent Water Molecules",
                                    "Total Persistent Hydrogen Bonds",
                                    "Total Persistent Hydrogen Bonds From Protein"
                                };
        d.fpp = xvgropen(opt2fn("-op", NFILE, fnm), "Persistent Waters and Hydrogen Bonds", "Number of Frames", "Number Water Molecules", oenv);
        xvgr_legend(d.fpp,asize(flegend), flegend, oenv);
    }
    
    /* Make sure -a1 and -a2 are included and increment by -1 to match internal numbering */
    if ( d.a1<0 || d.a2<0 ) {
        gmx_fatal(FARGS, "\nAtom numbers -a1 and -a2 defining the bond vector must be specified\n");
    }
    d.a1--; d.a2--;
    
    /* Make sure that -a1 and -a2 exist in the trajectory */
    if (d.a1<0 || d.a1>top->atoms.nr){
        gmx_fatal(FARGS, "\nError: Atom -a1 is %d, which is outside of the range 0 to %d\n",d.a1+1,top->atoms.nr+1);
    }
    if (d.a2<0 || d.a2>top->atoms.nr){
        gmx_fatal(FARGS, "\nError: Atom -a2 is %d, which is outside of the range 0 to %d\n",d.a2+1,top->atoms.nr+1);
    }

    /* Parse through the frames */
    gmx_ana_do(trj, 0, &analyze_frame, &d);
    
    /* Now we analyze the water information */
    analyze_information(&d) ; 

    if (d.bVerbose ) {
        fprintf(stderr, "Finished with analysis. Now closing all files\n") ; 
    }

    ffclose(d.fp); //This output file is automatically opened. Must close. 

    if (d.doLog) {
        ffclose(d.fpa);
    }
    if (d.doPersistent) {
        ffclose(d.fpp);
    }
    if (d.doWatGeo) {
        ffclose(d.fpg); 
    }
    if (d.doProtGeo) {
        ffclose(d.fpnwg); //Close all optional output files 
    }
    return 0 ; 
}


/*! \brief
 * The main function.
 *
 * In Gromacs, most analysis programs are implemented such that the \p main
 * function is only a wrapper for a \p gmx_something function that does all
 * the work, and that convention is also followed here. */
int main(int argc, char *argv[])
{
    gmx_nitrile_hbond(argc, argv);
    return 0;
}
