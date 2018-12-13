/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * scuff-scatter.cc  -- a standalone code within the scuff-EM suite for 
 *                   -- solving problems involving the scattering of 
 *                   -- electromagnetic radiation from an arbitrary 
 *                   -- compact object
 *
 * homer reid        -- 10/2006 -- 1/2012
 *
 * documentation at: http://homerreid.com/scuff-em/
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include "scuff-scatter.h"

#include <iostream>
/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXPS    10    // max number of point sources
#define MAXFREQ  10    // max number of frequencies
#define MAXEPF   10    // max number of evaluation-point files
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *scuff_scatter(char *GeoFile, char *OmegaFile, double *psLoc, int npsLoc, 
                    cdouble *psStrength, int npsStrength, char *EPFile, cdouble EH[6])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  InstallHRSignalHandler();
  InitializeLog((char*)"scuff-scatter");
  char *IFFile=0;
  char *TransFile=0;                 
  char *FileBase=0;

  HMatrix *PSDMatrix=0;
  if (FileBase==0) 
   FileBase=vstrdup(GetFileBase(GeoFile));

  /*******************************************************************/
  /* process frequency-related options                               */
  /*******************************************************************/
  HVector *OmegaList=GetOmegaList(OmegaFile, 0, 0,0, 0, 0);
  //printf("number of frequencies: %d\n",OmegaList->N);
  //printf("frequency: %f\n",OmegaList->GetEntry(OmegaList->N-1));
  /*******************************************************************/
  /* process incident-field-related options to construct the data    */
  /* used to generate the incident field in our scattering problem.  */
  /*******************************************************************/
 if ( npsLoc!=npsStrength )
   ErrExit("numbers of --psLocation and --psStrength options must agree");

  IncField *IF=0;
  for(int nps=0; nps<npsLoc; nps++)
   { PointSource *PS=new PointSource(psLoc + 3*nps, psStrength + 3*nps);
     PS->Next=IF;
     IF=PS;
   };

  IncFieldList *IFList=0;
  if (IFList!=0 && IF!=0)
   ErrExit("--IFFile is incompatible with other incident-field specifications");
  else if (IFFile!=0 && IF==0)
   IFList = ReadIncFieldList(IFFile);
  else if (IFFile==0 && IF!=0)
   IFList = AddIncFieldToList(IF,const_cast<char *>("Default"));


  /*******************************************************************/
  /* create the SSData structure containing everything we need to    */
  /* execute scattering calculations                                 */
  /*******************************************************************/
  SSData MySSData, *SSD=&MySSData;

  RWGGeometry *G      = SSD->G   = new RWGGeometry(GeoFile);
  HMatrix *M          = SSD->M   = G->AllocateBEMMatrix();
  HVector *RHS        = SSD->RHS = G->AllocateRHSVector();
  HVector *KN         = SSD->KN  = G->AllocateRHSVector();
  double *kBloch      = SSD->kBloch = 0;
  SSD->IF             = 0;
  SSD->TransformLabel = 0;
  SSD->IFLabel        = 0;
  SSD->FileBase       = FileBase;
  
  /*--------------------------------------------------------------*/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /*--------------------------------------------------------------*/
  GTCList GTCs=ReadTransFile(TransFile);
  int NumTransformations=GTCs.size();
  char *ErrMsg=G->CheckGTCList(GTCs);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  /*******************************************************************/
  /* for periodic geometries, all incident field sources that are    */
  /* active at a given time must involve  single incident field      */
  /* source, which must be a plane wave, and the bloch wavevector    */
  /* is extracted from the plane wave direction                      */
  /*******************************************************************/
  double kBlochBuffer[3];
  if (G->LDim>0)
   { if ( npsLoc!=0 )
      ErrExit("for extended geometries, the incident field must be a single plane wave");
     kBloch = SSD->kBloch = kBlochBuffer;
   };


  /*******************************************************************/
  /* if we have more than one geometrical transformation,            */
  /* allocate storage for BEM matrix blocks                          */
  /*******************************************************************/
  HMatrix **TBlocks=0, **UBlocks=0;
  int NS=G->NumSurfaces;
  if (NumTransformations>1)
   { int NADB = NS*(NS-1)/2; // number of above-diagonal blocks
     TBlocks  = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
     UBlocks  = (HMatrix **)mallocEC(NADB*sizeof(HMatrix *));
     for(int ns=0, nb=0; ns<NS; ns++)
      { 
        int nsMate = G->Mate[ns];
        if ( nsMate!=-1 )
         TBlocks[ns] = TBlocks[nsMate];
        else
         { int NBF=G->Surfaces[ns]->NumBFs;
           TBlocks[ns] = new HMatrix(NBF, NBF, M->RealComplex);
         };

        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int NBF=G->Surfaces[ns]->NumBFs;
           int NBFp=G->Surfaces[nsp]->NumBFs;
           UBlocks[nb] = new HMatrix(NBF, NBFp, M->RealComplex);
         };
      };
   };

  /*******************************************************************/
  /* loop over frequencies *******************************************/
  /*******************************************************************/
  
  for(int nFreq=0; nFreq<OmegaList->N; nFreq++)
   { 
     cdouble Omega = OmegaList->GetEntry(nFreq);
     SSD->Omega    = Omega;

     char OmegaStr[MAXSTR];
     z2s(Omega, OmegaStr);
     Log(" Working at frequency %s...",OmegaStr);

     /*******************************************************************/
     /* if we have more than one transformation, pre-assemble diagonal  */
     /* matrix blocks at this frequency; otherwise just assemble the    */
     /* whole matrix                                                    */
     /*******************************************************************/
     if (NumTransformations==1)
      G->AssembleBEMMatrix(Omega, kBloch, M);
     else
      for(int ns=0; ns<G->NumSurfaces; ns++)
       if (G->Mate[ns]==-1)
        G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, TBlocks[ns]);


     /*******************************************************************/
     /*******************************************************************/
     /*******************************************************************/
     for(int nt=0; nt<NumTransformations; nt++)
      {
        char TransformStr[100]="";
        if (TransFile)
         { G->Transform(GTCs[nt]);
           SSD->TransformLabel=GTCs[nt]->Tag;
           Log("Working at transformation %s...",SSD->TransformLabel);
           snprintf(TransformStr,100,"_%s",SSD->TransformLabel);
         };

        /*******************************************************************/
        /* assemble and insert off-diagonal blocks as necessary ************/
        /*******************************************************************/
        if (NumTransformations>1)
         { for(int ns=0, nb=0; ns<G->NumSurfaces; ns++)
            for(int nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
             G->AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch, UBlocks[nb]);

           for(int ns=0, nb=0; ns<G->NumSurfaces; ns++)
            { int RowOffset=G->BFIndexOffset[ns];
              M->InsertBlock(TBlocks[ns], RowOffset, RowOffset);
              for(int nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
               { int ColOffset=G->BFIndexOffset[nsp];
                 M->InsertBlock(UBlocks[nb], RowOffset, ColOffset);
                 M->InsertBlockAdjoint(UBlocks[nb], ColOffset, RowOffset);
               };
            };
         };

        /*******************************************************************/
        /* LU-factorize the BEM matrix to prepare for solving scattering   */
        /* problems                                                        */
        /*******************************************************************/
        Log("  LU-factorizing BEM matrix...");
        M->LUFactorize();

        /***************************************************************/
        /* loop over incident fields                                   */
        /***************************************************************/
        for(int nIF=0; nIF<IFList->NumIFs; nIF++)
         { 
           IF = SSD->IF = IFList->IFs[nIF];
           SSD->IFLabel = IFFile ? IFList->Labels[nIF] : 0;
           if (SSD->IFLabel)
            Log("  Processing incident field %s...",SSD->IFLabel);

           char IFStr[100]="";
           if (SSD->IFLabel)
            snprintf(IFStr,100,"_%s",SSD->IFLabel);
   
           /***************************************************************/
           /* assemble RHS vector and solve BEM system*********************/
           /***************************************************************/
           Log("  Assembling RHS vector...");
           G->AssembleRHSVector(Omega, kBloch, IF, KN);
           RHS->Copy(KN); // copy RHS vector for later 
           Log("  Solving the BEM system...");
           M->LUSolve(KN);
                    
           /***************************************************************/
           /* now process all requested outputs                           */
           /***************************************************************/
      
           /*--------------------------------------------------------------*/
           /*- panel source densities -------------------------------------*/
           /*--------------------------------------------------------------*/
           
           if (PSDMatrix==0)
             PSDMatrix=G->GetPanelSourceDensities(Omega, KN, 0); 
           else
             G->GetPanelSourceDensities(Omega, KN, PSDMatrix);
           
           /****************************************************************/
           /*--------------------------------------------------------------*/
           /*- scattered fields at user-specified points ------------------*/
           /*--------------------------------------------------------------*/
           
           HMatrix *XMatrix=new HMatrix(EPFile,1,"-ncol 3");
           HMatrix *SFMatrix = G->GetFields( 0, KN, Omega, kBloch, XMatrix); // scattered
           SFMatrix->GetEntries(0,":",EH);
           delete XMatrix;
         }; // for(int nIF=0; nIF<IFList->NumIFs; nIF++
      
        /*******************************************************************/
        /*******************************************************************/
        /*******************************************************************/
        G->UnTransform();

      }; // for(int nt=0; nt<NumTransformations; nt++)

   }; //  for(nFreq=0; nFreq<NumFreqs; nFreqs++)
  delete G;
  //printf("Thank you for your support.\n");
  return PSDMatrix; 
}
