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
typedef struct _seqlist
{
   char *sequence;
   int  seqLen;
   struct _seqlist *next;
}  SEQLIST;

   
   

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
REAL ScorePairwiseAlignmentIdentity(char *seq1, char *seq2);
REAL MeanWeightedDiversity(SEQLIST *seqlist, int nSeqs, int position);
REAL ScoreDiversity(char *seq1, char *seq2, int position);
REAL ScoreSimilarity(char *seq1, char *seq2, int position);
REAL MDMScore(char aa1, char aa2);



/************************************************************************/
int main(int argc, char **argv)
{
   return(0);
   
}

/************************************************************************/
REAL MeanWeightedDiversity(SEQLIST *seqlist, int nSeqs, int position)
{
   SEQLIST *s;
   int j, k;
   REAL Vr = (REAL)0.0,
      P_1j,
      v_1j;
   char *seq1, *seq2;
   

   seq1 = seqlist->sequence;
   
   for(j=1; j<nSeqs; j++)
   {
      s = seqlist;
      
      for(k=0; k<=j; k++)
      {
         NEXT(s);
      }
      seq2 = s->sequence;
      
      P_1j = ScorePairwiseAlignmentIdentity(seq1, seq2);
      v_1j = ScoreDiversity(seq1, seq2, position);
      
      Vr += (v_1j * P_1j);
   }

   Vr /= (nSeqs - 1);

   return(Vr);
}

/************************************************************************/
REAL ScoreDiversity(char *seq1, char *seq2, int position)
{
   REAL similarity;
   
   similarity = ScoreSimilarity(seq1, seq2, position);
   return(1-similarity);
}

/************************************************************************/
REAL ScoreSimilarity(char *seq1, char *seq2, int position)
{
   return(MDMScore(seq1[position], seq2[position])/10.0);
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



REAL MDMScore(char aa1, char aa2)
{
   return(0);
}
   
