/************************************************************************/
/**

   Program:    carcons
   \file       carcons.c
   
   \version    V0.1
   \date       08.02.23  
   \brief      Score Conservation with Allowed Resiudes
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2023
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               Division of Biosciences,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1    03.02.23  Original  By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include "bioplib/MathType.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
REAL ScorePairwiseAlignmentIdentity(char *seq1, char *seq2);


/************************************************************************/
int main(int argc, char **argv)
{
   return(0);
   
}

/************************************************************************/
REAL ScorePairwiseAlignmentIdentity(char *seq1, char *seq2)
{
   int NPos   = 0,
      NMatch = 0,
      seqLen,
      i;

   seqLen = MIN(strlen(seq1), strlen(seq2));

   for(i=0; i<seqLen; i++)
   {
      if(seq1[i] == seq2[i])
      {
         /* It's a match                                                */
         if(seq1[i] != '-')     /* It's not a double insertion          */
         {
            NPos++;
            NMatch++;
         }
      }
      else
      {
         /* It's not a match                                            */
         NPos++;
      }
   }
   
   return((REAL)NMatch / (REAL)NPos);
}
