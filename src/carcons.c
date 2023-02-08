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
   Compile with -std=c99 to activate inlining


**************************************************************************

   Revision History:
   =================
   V0.1    03.02.23  Original  By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "bioplib/sequtil.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 512

typedef struct _seqlist
{
   char *sequence;
   int  seqLen;
   struct _seqlist *next;
}  SEQLIST;

/* If compiled with -std=c99 some functions will be inlined             */
#ifdef __GNUC_STDC_INLINE__
#   define INLINE static inline
#else
#   define INLINE static
#endif

   

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
REAL ScorePairwiseAlignmentIdentity(char *seq1, char *seqJ);
REAL SDWeightedDiversity(SEQLIST *seqlist, int nSeqs, int position,
                         REAL meanWeightedDiversity);
REAL MeanWeightedDiversity(SEQLIST *seqlist, int nSeqs, int position);
INLINE REAL ScoreDiversity(char *seq1, char *seqJ, int position);
INLINE REAL ScoreSimilarity(char *seq1, char *seqJ, int position);
REAL MDMScore(char aa1, char aa2);
INLINE REAL CalculateThreshold(REAL meanWeightedDiversity,
                        REAL sdWeightedDiversity,
                        REAL numSD);
REAL ScoreSelfIDMutant(char *seq, int nMutant);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  int *position, char *mutation, REAL *numSD);
void Usage(void);
SEQLIST *ReadFASTAAlignment(FILE *in, int *nSeqs);
REAL EvaluateMutation(SEQLIST *seqlist, int NSeqs,
                      int position, char mutation, REAL numSD);
REAL CalculateDelta(SEQLIST *seqlist, int position, char mutant,
                    REAL threshold);





/************************************************************************/
int main(int argc, char **argv)
{
   FILE    *in  = stdin,
           *out = stdout;
   char    inFile[MAXBUFF],
           outFile[MAXBUFF],
           mutation;
   REAL    numSD    = 2.0;
   int     position = 0,
           nSeqs    = 0;
   SEQLIST *seqlist = NULL;
   
   if(ParseCmdLine(argc, argv, inFile, outFile, &position, &mutation,
                   &numSD))
   {
      if(blOpenStdFiles(inFile, outFile, &in, &out))
      {
         if((seqlist = ReadFASTAAlignment(in, &nSeqs))!=NULL)
         {
            REAL score = EvaluateMutation(seqlist, nSeqs,
                                          position, mutation, numSD);

            fprintf(out, "%.3f\n", score);
         }
         else
         {
            fprintf(stderr,"carcons: Error - unable to read FASTA \
alignment file\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"carcons: Error - unable to open input or \
output file\n");
         return(1);

      }
   }
   else
   {
      Usage();
   }
   
   return(0);
   
}

REAL EvaluateMutation(SEQLIST *seqlist, int nSeqs,
                      int position, char mutation, REAL numSD)
{
   REAL meanWeightedDiversity,
        sd,
        threshold,
        score;

   meanWeightedDiversity = MeanWeightedDiversity(seqlist,nSeqs,position);
   sd                    = SDWeightedDiversity(seqlist,nSeqs,position,
                                               meanWeightedDiversity);
   threshold             = CalculateThreshold(meanWeightedDiversity,
                                              sd, numSD);
   score                 = CalculateDelta(seqlist, position, mutation,
                                          threshold);

   return(score);
}


/************************************************************************/
SEQLIST *ReadFASTAAlignment(FILE *in, int *nSeqs)
{
   SEQLIST *seqlist = NULL,
           *s;
   char    *sequence,
           header[MAXBUFF];

   *nSeqs = 0;
   
   while((sequence = blReadFASTA(in, header, MAXBUFF))!=NULL)
   {
      if(seqlist == NULL)
      {
         INIT(seqlist, SEQLIST);
         s = seqlist;
      }
      else
      {
         ALLOCNEXT(s, SEQLIST);
      }

      if(s==NULL)
      {
         FREELIST(seqlist, SEQLIST);
         return(NULL);
      }

      s->sequence = sequence;
      s->seqLen   = strlen(sequence);
      (*nSeqs)++;
   }
   
   return(seqlist);
}


/************************************************************************/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  int *position, char *mutation, REAL *numSD)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   if(argc < 2)
      return(FALSE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'n':
            argc--;
            argv++;
            if(!sscanf(argv[0], "%lf", numSD))
               return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2-4 <= 2 arguments left                */
         if((argc < 2) || (argc > 4))
            return(FALSE);

         /* Copy the first (required) to position                       */
         if(!sscanf(argv[0], "%d", position))
            return(FALSE);
         argc--;
         argv++;

         /* Copy the second (required) to mutation                      */
         *mutation = argv[0][0];
        
         /* If there is another, copy the first to infile               */
         argc--;
         argv++;
         if(argc)
         {
            strcpy(infile, argv[0]);
         }
         
         /* If there is another, copy the second to outfile             */
         argc--;
         argv++;
         if(argc)
         {
            strcpy(infile, argv[0]);
         }
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}



/************************************************************************/
REAL MeanWeightedDiversity(SEQLIST *seqlist, int nSeqs, int position)
{
   SEQLIST *s;
   int     j, k;
   REAL    Vr = (REAL)0.0,
           P_1j,
           v_1j;
   char    *seq1, *seqJ;
   

   seq1 = seqlist->sequence;
   
   for(j=1; j<nSeqs; j++)
   {
      s = seqlist;
      
      for(k=0; k<j; k++)
      {
         NEXT(s);
      }
      seqJ = s->sequence;
      
      P_1j = ScorePairwiseAlignmentIdentity(seq1, seqJ);
      v_1j = ScoreDiversity(seq1, seqJ, position);
      
      Vr += (v_1j * P_1j);
   }

   Vr /= (nSeqs - 1);

   return(Vr);
}

/************************************************************************/
INLINE REAL CalculateThreshold(REAL meanWeightedDiversity,
                               REAL sdWeightedDiversity,
                               REAL numSD)
{
   return(meanWeightedDiversity + (numSD * sdWeightedDiversity));
}



/************************************************************************/
REAL CalculateDelta(SEQLIST *seqlist, int position, char mutant,
                    REAL threshold)
{
   REAL P_1M,
        v_1M;
   char native;

   native = seqlist->sequence[position];

   v_1M = 1 - (MDMScore(native, mutant)/10.0);
   P_1M = ScoreSelfIDMutant(seqlist->sequence, 1);
   
   return(threshold - (v_1M * P_1M));
}


/************************************************************************/
REAL SDWeightedDiversity(SEQLIST *seqlist, int nSeqs, int position,
                         REAL meanWeightedDiversity)
{
   char *seq1, *seqJ;
   int  j, k;
   
   REAL sumDevSq = 0.0,
        dev, sd,
        P_1j,
        v_1j;
   seq1 = seqlist->sequence;
   
   for(j=1; j<nSeqs; j++)
   {
      SEQLIST *s;
      
      s = seqlist;
      for(k=0; k<j; k++)
      {
         NEXT(s);
      }
      seqJ = s->sequence;
      
      P_1j = ScorePairwiseAlignmentIdentity(seq1, seqJ);
      v_1j = ScoreDiversity(seq1, seqJ, position);

      dev = ((v_1j * P_1j) - meanWeightedDiversity);
      sumDevSq += (dev*dev);
   }

   sd = sqrt(sumDevSq/(nSeqs-1));

   return(sd);
}



/************************************************************************/
INLINE REAL ScoreDiversity(char *seq1, char *seqJ, int position)
{
   REAL similarity;
   
   similarity = ScoreSimilarity(seq1, seqJ, position);
   return(1-similarity);
}

/************************************************************************/
INLINE REAL ScoreSimilarity(char *seq1, char *seqJ, int position)
{
   return(MDMScore(seq1[position], seqJ[position])/10.0);
}

/************************************************************************/
REAL ScorePairwiseAlignmentIdentity(char *seq1, char *seqJ)
{
   int nPos   = 0,
       nMatch = 0,
       seqLen,
       i;

   seqLen = MIN(strlen(seq1), strlen(seqJ));

   for(i=0; i<seqLen; i++)
   {
      if(seq1[i] == seqJ[i])
      {
         /* It's a match                                                */
         if(seq1[i] != '-')     /* It's not a double insertion          */
         {
            nPos++;
            nMatch++;
         }
      }
      else
      {
         /* It's not a match                                            */
         nPos++;
      }
   }
   
   return((REAL)nMatch / (REAL)nPos);
}
/************************************************************************/
REAL ScoreSelfIDMutant(char *seq, int nMutant)
{
   int nPos   = 0,
       nMatch = 0,
       seqLen,
       i;

   seqLen = strlen(seq);

   for(i=0; i<seqLen; i++)
   {
      if(seq[i] != '-')     /* It's not an insertion                    */
      {
         nPos++;
         nMatch++;
      }
   }
   
   return((REAL)(nMatch-nMutant) / (REAL)nPos);
}

void Usage(void)
{
   printf("Usage\n");
}



REAL MDMScore(char aa1, char aa2)
{
   return(0);
}
   
