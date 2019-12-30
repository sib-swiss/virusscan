/*	------------------------------------------------------------------------------------

                                  * virusscan *
	Detection and quantification of virus sequences from illumina fastq files


    Copyright (C) SIB  - Swiss Institute of Bioinformatics,   2015-2019 Nicolas Guex and Christian Iseli
    Copyright (C) UNIL - University of Lausanne, Switzerland       2019 Nicolas Guex and Christian Iseli


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


	Code:       Nicolas Guex and Christian Iseli
	Contacts:   Nicolas.Guex@unil.ch and Christian.Iseli@unil.ch
	Repository: https://github.com/sib-swiss/virusscan



	Machine :	Unix
	Language:	C
	Requires:	AVX512, pthread

	Version information

	Version:	1.0  Dec.  2019 Public release of code under GPL2+ license

	
	References:  
	
	      http://clinicalmetagenomics.org/wp-content/uploads/2017/11/Lemercier.pdf

	     to obtain the virus database compatible with this release (20180906.bin)
		 contact P. Lemercier (philippe.lemercier@sib.swiss)
		 
		 a newer version, based on the latest Virosaurus database will come in 2020
		               https://viralzone.expasy.org/8676

 

	Compiling:

	gcc  -mavx512f  -O3  -o virus_match_AVX512 src/virus_match.c  -lpthread

	------------------------------------------------------------------------------------
*/

//=================================================================================================

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <semaphore.h>
#include <pthread.h>
#include <assert.h>
#include <immintrin.h>
#include <errno.h>
#include <termios.h>
#include <sys/ioctl.h>

#include <sys/io.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <signal.h>

#include <sys/ipc.h>
#include <sys/sem.h>


#include "virus.h"



#define k128bitsLongAccess 16L

//=================================================================================================
#define kMaxFileName 512
#define kTagSizeStorage 256
#define kMaxInputLineLen 1024  
#define kMAX_TAG_SNDRCV_THREAD_cnt 512  
#define kINPUTBUFFER_MAXCNT_PER_CPU 128L
#define kINPUTBUFFER_MAXCNT (kMAX_TAG_SNDRCV_THREAD_cnt * kINPUTBUFFER_MAXCNT_PER_CPU)
#define kINPUTBUFFER_DATA_SIZE  (kINPUTBUFFER_MAXCNT*kMaxInputLineLen)
#define kMaxMultGeneHit 512 
#define kMaxHitsPerCPU 32768

#define kMaxMismatchPerTagDBcheck 15
#define kMaxNperTag 50
#define kMaxMismatchPerTag 12
#define kMaxMismatchPerPair 15 

//=================================================================================================

#define kFWD 0
#define kREV 1

//---------------------------------------------------------------

typedef	struct MAPPINGCNT_struct  MAPPINGCNT;
struct MAPPINGCNT_struct
{
	unsigned int cnt;
	unsigned int maxcnt;
};

//---------------------------------------------------------------

typedef	struct VIROSAURUSRSLT_struct  VIROSAURUSRSLT;
struct VIROSAURUSRSLT_struct
{
	unsigned long virosauruspos;
	unsigned int tagid;
	unsigned int tagid2;
	unsigned int taglen;
	         int shift;
	unsigned int tag2shift;
	unsigned int tag2len;
	         int flags;
	unsigned int mismatch;
	unsigned int splitAssignment;  
	unsigned int multipleVirusHit;  
	unsigned int isMappingOnSingleViruskind;
	unsigned int UID;
};

//---------------------------------------------------------------

typedef	struct THREAD_DATA_struct  THREAD_DATA;
struct THREAD_DATA_struct
{
	unsigned int cpu;
	unsigned int filterVirusID;
};

//---------------------------------------------------------------
typedef	struct VIRUSCNT_struct  VIRUSCNT;
struct VIRUSCNT_struct
{
	float virus;
	float virusUniqueGood;
	unsigned int  nameidx;
};

//---------------------------------------------------------------

static 	unsigned char *volatile inputbuffer;
static  unsigned long volatile gProcessedCnt[kMAX_TAG_SNDRCV_THREAD_cnt];
static  unsigned long volatile gCrappy[kMAX_TAG_SNDRCV_THREAD_cnt];
static  unsigned long volatile gUnmappedCnt[kMAX_TAG_SNDRCV_THREAD_cnt];


static struct timeval ts;
static struct timeval te;

static unsigned char gNT1bit[128];
static unsigned char grevNT1bit[128];

static 	unsigned int *indexLeft = NULL;
static char *tagNames1 = NULL;
static unsigned long *tagUID1 = NULL;
static unsigned long *virosaurusDBidx = NULL;
static char *virosaurusDB = NULL;
static	unsigned char *sortedtagsLeft = NULL;
static	HASH *hashLeft = NULL;
static char *virusshortnames = NULL;

static char gIUPACrevcomp[128];


static char gOutputDir[kMaxFileName];
static char gInputFile[kMaxFileName];
static char gInputFile2[kMaxFileName];
static unsigned long gMaxUID = 0L;

static int gThreadCnt = 1;
static unsigned int gSingleEnd = 0;
static FILE *rsltfile[kMAX_TAG_SNDRCV_THREAD_cnt];

#define kMAPPINGRESULTS_CHUNK 32768
static MAPPINGCNT gMappingResultsCnt[kMAX_TAG_SNDRCV_THREAD_cnt][kMAX_VIRUS_CHUNKS];
static VIROSAURUSRSLT *gMappingResults[kMAX_TAG_SNDRCV_THREAD_cnt][kMAX_VIRUS_CHUNKS];


static  char gFilterAdapter[256];
static VIRUSCNT *gVirusCnt;

static unsigned int gMinTagLen;
static unsigned int gMinFragLen;
static unsigned int gMaxMismatchForUnambiguous;
static int gMostFreqFragLen = 418;
//=================================================================================================
//=================================================================================================
//=================================================================================================

void initRevCompTable(unsigned char *gIUPACrevcomp)
{

        memset(gIUPACrevcomp,'N',128);
        gIUPACrevcomp['A'] = 'T';
        gIUPACrevcomp['T'] = 'A';
        gIUPACrevcomp['U'] = 'A';
        gIUPACrevcomp['G'] = 'C';
        gIUPACrevcomp['C'] = 'G';
        gIUPACrevcomp['Y'] = 'R';
        gIUPACrevcomp['R'] = 'Y';
        gIUPACrevcomp['S'] = 'S';
        gIUPACrevcomp['W'] = 'W';
        gIUPACrevcomp['K'] = 'M';
        gIUPACrevcomp['M'] = 'K';
        gIUPACrevcomp['B'] = 'V';
        gIUPACrevcomp['D'] = 'H';
        gIUPACrevcomp['H'] = 'D';
        gIUPACrevcomp['V'] = 'B';

        gIUPACrevcomp['a'] = 't';
        gIUPACrevcomp['t'] = 'a';
        gIUPACrevcomp['u'] = 'a';
        gIUPACrevcomp['g'] = 'c';
        gIUPACrevcomp['c'] = 'g';
        gIUPACrevcomp['y'] = 'r';
        gIUPACrevcomp['r'] = 'y';
        gIUPACrevcomp['s'] = 's';
        gIUPACrevcomp['w'] = 'w';
        gIUPACrevcomp['k'] = 'm';
        gIUPACrevcomp['m'] = 'k';
        gIUPACrevcomp['b'] = 'v';
        gIUPACrevcomp['d'] = 'h';
        gIUPACrevcomp['h'] = 'd';
        gIUPACrevcomp['v'] = 'b';
        gIUPACrevcomp['n'] = 'n';

} // initRevCompTable
//---------------------------------------------------------------

static unsigned int packFwdTag(char *tp,unsigned char *packed4nt,unsigned int l)
{
	unsigned char *packed4ntPtr;
	unsigned int pos;
	unsigned int shift;
	unsigned int encoded = 1;
		packed4ntPtr = &packed4nt[0];

	__m512i	ZEROvec = _mm512_xor_epi32(ZEROvec,ZEROvec);
	_mm512_store_epi32(packed4ntPtr,ZEROvec);


		pos = 0;
		shift = 0;
		pos = 0;
		while(pos < l)
		{
			*packed4ntPtr |= (gNT1bit[*tp] << shift);
			shift +=1;
			if (shift == 8)
			{
				shift = 0;
				packed4ntPtr++;
			}
			tp++;
			pos++;
		}

		return(encoded);
	
} // packFwdTag
//---------------------------------------------------------------

static int openPackedTags(char *fn, HEADER *hdr)
{
	FILE *f = NULL;

		f = fopen(fn,"rb");
		if (f == NULL)
		{
			printf("Error reading %s\n",fn);
			return(-1);
		}
		fread(hdr,sizeof(HEADER),1,f);
		fclose(f);
		
		return(0);
	
} // openPackedTags
//---------------------------------------------------------------

static int loadVirosaurus(unsigned long *tagUID,char *tagNames, unsigned char *sortedtagsLeft,unsigned int *indexLeft,HASH *hashLeft,unsigned int tagcnt,char *virosaurusDB,unsigned long *virosaurusDBidx,char *fn)
{
	FILE *f = NULL;
	HEADER hdr;	

		f = fopen(fn,"rb");
		if (f == NULL)
		{
			printf("Error reading %s\n",fn); fflush(stdout);
			return(-1);
		}

		fread(&hdr,sizeof(HEADER),1,f);
		
		if (hdr.version != 11)
		{
			printf("Wrong version format %s (v%d but expecting v11)\n",fn,hdr.version); fflush(stdout);
			return(-1);
		}
		if (hdr.hashSize != kHashSize)
		{
			printf("Wrong hash format %s\n",fn); fflush(stdout);
			return(-1);
		}
		if (hdr.taglen != kDatabaseTagStorageSize)
		{
			printf("Wrong tag size %s %d but compiled for %d\n",fn,hdr.taglen,kDatabaseTagStorageSize); fflush(stdout);
			return(-1);
		}
		
		
		fread(sortedtagsLeft,k128bits*sizeof(unsigned char),tagcnt,f);
		fread(indexLeft,sizeof(unsigned int),tagcnt,f);
		fread(hashLeft,sizeof(HASH),kHashSize,f);

		fread(&tagNames[0           ],kTagNameMaxLen*sizeof(unsigned char),tagcnt,f);
		fread(&tagUID[0           ],sizeof(unsigned long),tagcnt,f);

		
		fread(&virosaurusDBidx[0           ],sizeof(unsigned long),tagcnt,f);
		fread(&virosaurusDB[0           ],sizeof(char),hdr.virosaurusDBsize,f);

		fread(&virusshortnames[0    ],kShortNameLen*sizeof(unsigned char),hdr.maxVirusID+1,f);

		fclose(f);
	
		return(0);
	
} // loadVirosaurus
//---------------------------------------------------------------

static unsigned char gNT1bit[128];
static unsigned char grevNT1bit[128];
static unsigned char *packedtags = NULL;
static HASH	*memhash = NULL;
static unsigned int *mem512 = NULL;
static unsigned int *masks512 = NULL;


static void initPackNT1bitTable(unsigned char *gNT1bit)
{
	memset(gNT1bit,0,128);
	gNT1bit['A'] = 0x0;
	gNT1bit['C'] = 0x0;
	gNT1bit['G'] = 0x1;
	gNT1bit['T'] = 0x1;

} // initPackNT1bitTable
//---------------------------------------------------------------
static void initPackRevNT1bitTable(unsigned char *grevNT1bit)
{
	memset(grevNT1bit,0x1 ,128);
	grevNT1bit['A'] = 0x1;
	grevNT1bit['C'] = 0x1;
	grevNT1bit['G'] = 0x0;
	grevNT1bit['T'] = 0x0;

} // initPackRevNT1bitTable
//---------------------------------------------------------------

#define kCopyLeftToAll 0
static void init512bitsVectors(void)
{
	unsigned int i;
	unsigned int offset;
	
	offset = kCopyLeftToAll;
	for (i = 0; i< 4; i++)
	{
		mem512[offset+i] = i;
		mem512[offset+i+4] = i;
		mem512[offset+i+8] = i;
		mem512[offset+i+12] = i;
	}

} // init512bitsVectors
//---------------------------------------------------------------


static void init512bitsMasks(void)
{
	unsigned int i,l;
	unsigned int offset = 0;
	char tagMask[kTagSizeStorage];
	__m512i tmp;
	__m512i cpy = _mm512_load_epi32(&mem512[kCopyLeftToAll]);
	
	for (l = 0; l < kTagSizeStorage; l++)
	{

		
		for (i = 0; i < l; i++)
			tagMask[i] = 'G';				
		for (i = l; i < kDatabaseTagStorageSize; i++)
			tagMask[i] = 'A';				


		packFwdTag(tagMask ,(unsigned char *)&masks512[offset] ,l);
		tmp = _mm512_load_epi32(&masks512[offset]);
		tmp = _mm512_permutexvar_epi32(cpy,tmp);
		_mm512_store_epi32(&masks512[offset],tmp);
		offset += 16;  

	}

} // init512bitsMasks

//---------------------------------------------------------------

static unsigned int packVarLenFwdTag(char *tp,unsigned char *packed4nt)
{
	unsigned char *packed4ntPtr;
	unsigned int pos;
	unsigned int shift;

		packed4ntPtr = &packed4nt[0];

		__m512i	ZEROvec = _mm512_xor_epi32(ZEROvec,ZEROvec);
		_mm512_store_epi32(packed4ntPtr,ZEROvec);

		pos = 0;
		shift = 0;
		pos = 0;


		while(*tp != 0)
		{

			*packed4ntPtr |= (gNT1bit[*tp] << shift);
			shift +=1;
			if (shift == 8)
			{
				shift = 0;
				packed4ntPtr++;
			}
			tp++;
			pos++;
			if (pos >= kDatabaseTagStorageSize)
			{
				while(*tp++ != 0) { pos++; };
				break;
			}
		}

		return(pos);
	
} // packVarLenFwdTag
//---------------------------------------------------------------

static void packVarLenRevTag(char *tp,char *revtag, unsigned char *packed4nt,unsigned int l)
{
	unsigned char *packed4ntPtr;
	unsigned int pos;
	unsigned int shift;

		packed4ntPtr = &packed4nt[0];

		__m512i	ZEROvec = _mm512_xor_epi32(ZEROvec,ZEROvec);
		_mm512_store_epi32(packed4ntPtr,ZEROvec);

		pos = 0;
		shift = 0;
		pos = 0;

		while(pos < l)
		{
			revtag[pos] = gIUPACrevcomp[*tp];

			*packed4ntPtr |= (grevNT1bit[*tp] << shift);
			shift +=1;
			if (shift == 8)
			{
				shift = 0;
				packed4ntPtr++;
			}
			tp--;
			pos++;
			if (pos >= kDatabaseTagStorageSize) 
				break;
		}

		revtag[pos] = 0;
	
} // packVarLenRevTag
//---------------------------------------------------------------

unsigned int DoSearchVirosaurus2(char *tag1,char *revtag1,unsigned char *queryTag512bits,VIROSAURUSRSLT *rslts,unsigned int *minMM)
{
		unsigned int hash;
		unsigned int starttagpairid;
		unsigned int laststarttagpairid;
		unsigned int taglen;
		unsigned int i;
		unsigned int hitcnt = 0;
		int shift = 0;
		unsigned int minMismatch = 256;

		__int64_t bitCount[8] __attribute__((aligned(64)));

		__m512i A1vec;  
		__m512i B1vec;	
		__m512i vec;	
		__m512i ZEROvec;
		__m512i tvec;
		__m512i shiftmask;
		__m512i cpy = _mm512_load_epi32(&mem512[kCopyLeftToAll]);

		ZEROvec = _mm512_xor_epi32(ZEROvec,ZEROvec);

		
		for (shift = 36; shift >= 0; shift -=18)
		{
			taglen = packVarLenFwdTag(&tag1[shift]           ,  &queryTag512bits[0]);
			tvec = _mm512_load_epi32(&queryTag512bits[0]);

			shiftmask = _mm512_load_epi32(&masks512[taglen*16]);   
			A1vec = _mm512_permutexvar_epi32 (cpy,tvec);
			hash = *((unsigned int *)&queryTag512bits[kHashEncodingPos]);
			hash >>= kHashShift;
			hash &= kHashMask;


			starttagpairid = memhash[hash].fromPos ;
			starttagpairid &= 0xFFFFFFFC; 

			laststarttagpairid = memhash[hash].toPos+4;


			for (i = starttagpairid;i <= laststarttagpairid; i+=4)
			{
				
				B1vec = _mm512_load_epi32(&packedtags[i*k128bitsLongAccess]);

				
				vec = _mm512_xor_epi32(A1vec,B1vec);
				vec = _mm512_and_epi32(vec,shiftmask);
				
				_mm512_store_epi32(&bitCount[0], vec);
				unsigned int mmCnt;
				mmCnt = _mm_popcnt_u64(bitCount[0]) + _mm_popcnt_u64(bitCount[1]);
				if (mmCnt <= kMaxMismatchPerTagDBcheck)
				{
					rslts[hitcnt].tagid = i;
					rslts[hitcnt].taglen = taglen;
					rslts[hitcnt].shift = shift;
					rslts[hitcnt].mismatch = mmCnt;
					rslts[hitcnt++].flags = kFWD;
					if (mmCnt < minMismatch)
						minMismatch = mmCnt;
				}
				mmCnt = _mm_popcnt_u64(bitCount[2]) + _mm_popcnt_u64(bitCount[3]);
				if (mmCnt <= kMaxMismatchPerTagDBcheck)
				{
					rslts[hitcnt].tagid = i+1;
					rslts[hitcnt].taglen = taglen;
					rslts[hitcnt].shift = shift;
					rslts[hitcnt].mismatch = mmCnt;
					rslts[hitcnt++].flags = kFWD;
					if (mmCnt < minMismatch)
						minMismatch = mmCnt;
				}
				mmCnt = _mm_popcnt_u64(bitCount[4]) + _mm_popcnt_u64(bitCount[5]);
				if (mmCnt <= kMaxMismatchPerTagDBcheck)
				{
					rslts[hitcnt].tagid = i+2;
					rslts[hitcnt].taglen = taglen;
					rslts[hitcnt].shift = shift;
					rslts[hitcnt].mismatch = mmCnt;
					rslts[hitcnt++].flags = kFWD;
					if (mmCnt < minMismatch)
						minMismatch = mmCnt;
				}
				mmCnt = _mm_popcnt_u64(bitCount[6]) + _mm_popcnt_u64(bitCount[7]);
				if (mmCnt <= kMaxMismatchPerTagDBcheck)
				{
					rslts[hitcnt].tagid = i+3;
					rslts[hitcnt].taglen = taglen;
					rslts[hitcnt].shift = shift;
					rslts[hitcnt].mismatch = mmCnt;
					rslts[hitcnt++].flags = kFWD;
					if (mmCnt < minMismatch)
						minMismatch = mmCnt;
				}
				if (hitcnt >= (kMaxHitsPerCPU-16))
					break;
				
			} 

		
			packVarLenRevTag(&tag1[taglen-1],&revtag1[0],&queryTag512bits[0]  ,taglen);
			tvec = _mm512_load_epi32(&queryTag512bits[0]); 
			shiftmask = _mm512_load_epi32(&masks512[taglen*16]);   

			

			A1vec = _mm512_permutexvar_epi32 (cpy,tvec);

			hash = *((unsigned int *)&queryTag512bits[kHashEncodingPos]);  
			hash >>= kHashShift;
			hash &= kHashMask;

			starttagpairid = memhash[hash].fromPos ;
			starttagpairid &= 0xFFFFFFFC; 

			laststarttagpairid = memhash[hash].toPos+4;


			for (i = starttagpairid;i <= laststarttagpairid; i+=4)
			{
				
				B1vec = _mm512_load_epi32(&packedtags[i*k128bitsLongAccess]);

				
				vec = _mm512_xor_epi32(A1vec,B1vec);
				vec = _mm512_and_epi32(vec,shiftmask);
				_mm512_store_epi32(&bitCount[0], vec);
				unsigned int mmCnt;
				mmCnt = _mm_popcnt_u64(bitCount[0]) + _mm_popcnt_u64(bitCount[1]);
				if (mmCnt <= kMaxMismatchPerTagDBcheck)
				{
					rslts[hitcnt].tagid = i;
					rslts[hitcnt].taglen = taglen;
					rslts[hitcnt].shift = shift;
					rslts[hitcnt].mismatch = mmCnt;
					rslts[hitcnt++].flags = kREV;
					if (mmCnt < minMismatch)
						minMismatch = mmCnt;
				}
				mmCnt = _mm_popcnt_u64(bitCount[2]) + _mm_popcnt_u64(bitCount[3]);
				if (mmCnt <= kMaxMismatchPerTagDBcheck)
				{
					rslts[hitcnt].tagid = i+1;
					rslts[hitcnt].taglen = taglen;
					rslts[hitcnt].shift = shift;
					rslts[hitcnt].mismatch = mmCnt;
					rslts[hitcnt++].flags = kREV;
					if (mmCnt < minMismatch)
						minMismatch = mmCnt;
				}
				mmCnt = _mm_popcnt_u64(bitCount[4]) + _mm_popcnt_u64(bitCount[5]);
				if (mmCnt <= kMaxMismatchPerTagDBcheck)
				{
					rslts[hitcnt].tagid = i+2;
					rslts[hitcnt].taglen = taglen;
					rslts[hitcnt].shift = shift;
					rslts[hitcnt].mismatch = mmCnt;
					rslts[hitcnt++].flags = kREV;
					if (mmCnt < minMismatch)
						minMismatch = mmCnt;
				}
				mmCnt = _mm_popcnt_u64(bitCount[6]) + _mm_popcnt_u64(bitCount[7]);
				if (mmCnt <= kMaxMismatchPerTagDBcheck)
				{
					rslts[hitcnt].tagid = i+3;
					rslts[hitcnt].taglen = taglen;
					rslts[hitcnt].shift = shift;
					rslts[hitcnt].mismatch = mmCnt;
					rslts[hitcnt++].flags = kREV;
					if (mmCnt < minMismatch)
						minMismatch = mmCnt;
				}

				if (hitcnt >= (kMaxHitsPerCPU-16))
					break;
				
			} 

		}
		*minMM = minMismatch;
		return(hitcnt);
	
} // DoSearchVirosaurus2
//---------------------------------------------------------------

unsigned int VirosaurusCheckHits2(VIROSAURUSRSLT *rslts, unsigned int rsltcnt, unsigned int *index, char *fwd, char *rev,unsigned int minMismatch)
{
		int r;
		unsigned int rok = 0;
		unsigned int k;
	
	
		for (r = 0; r < rsltcnt; r++)
		{
			if (rslts[r].mismatch > (minMismatch+3))
				continue;


			unsigned long tp = index[ rslts[r].tagid ];
			if (virosaurusDBidx[tp] <  (unsigned long)rslts[r].shift)
			{
				goto alreadyseen;
			}
			

			for (k = 0; k < rok; k++)
			{
				if ((virosaurusDBidx[tp] - rslts[r].shift) == (virosaurusDBidx[rslts[k].tagid] - rslts[k].shift))
				{
				
					goto alreadyseen;
				}
			}

			
			unsigned int mismatch = 0;
			unsigned int mi = 0;
			unsigned long tidx = virosaurusDBidx[tp] - rslts[r].shift;
			unsigned int ti = 0;
			unsigned int Ncnt = 0;
			if (rslts[r].flags & kREV)
			{
				while(rev[mi] != 0)
				{
					if (rev[mi] != virosaurusDB[tidx+ti])
					{
						if ((rev[mi] == 'A') || (rev[mi] == 'C') || (rev[mi] == 'G') || (rev[mi] == 'T'))
							mismatch++;
						else
							Ncnt++;
						
						
					}
					else if (rev[mi] == 'N')
						Ncnt++;

					mi++;
					ti++;
				}
			}
			else
			{
				while(fwd[mi] != 0)
				{
					if (fwd[mi] != virosaurusDB[tidx+ti])
					{
						if ((fwd[mi] == 'A') || (fwd[mi] == 'C') || (fwd[mi] == 'G') || (fwd[mi] == 'T'))
							mismatch++;
						else
							Ncnt++;
						
						
					}
					else if (fwd[mi] == 'N')
						Ncnt++;
					
					mi++;
					ti++;
				}
			}
		
		

			if (mismatch <= kMaxMismatchPerTag && Ncnt <= kMaxNperTag)
			{
				rslts[rok].tagid = tp;
				rslts[rok].mismatch = mismatch;
				rslts[rok].virosauruspos = tagUID1[tp];
				rslts[rok].flags = rslts[r].flags;
				rslts[rok].shift = rslts[r].shift;
				rslts[rok].tag2len = mi;
				rok++;
			}
			alreadyseen:;
		}
		return(rok);
	
} // VirosaurusCheckHits2
//---------------------------------------------------------------
int sortMultVirusFunc(const void * a, const void * b)
{
   return ( *(int*)b - *(int*)a );
}
//---------------------------------------------------------------
int ValidatePairs(FILE *f,unsigned long UID,VIROSAURUSRSLT *rslts1, unsigned int rsltcnt1,VIROSAURUSRSLT *rslts2, unsigned int rsltcnt2,char *tag1,char *tag2,char *revtag1,char *revtag2, unsigned int cpu)
{
	unsigned int i1,i2;
	VIROSAURUSRSLT	validrslts[kMaxHitsPerCPU];
	unsigned int validrsltsCnt = 0;
	unsigned int bestmismatch = 999999;

	
	// first pass to determine lowest mismatch
	for (i1 = 0; i1 < rsltcnt1; i1++)
	{
		unsigned int virus1 = (unsigned int)((rslts1[i1].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT));
		unsigned int offset1 = (unsigned int)(rslts1[i1].virosauruspos & kPOSITION_MASK);
		for (i2 = 0; i2 < rsltcnt2; i2++)
		{
			unsigned int virus2 = (unsigned int)((rslts2[i2].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT));
			if (virus1 == virus2)	
			{

				unsigned int offset2 = (unsigned int)(rslts2[i2].virosauruspos & kPOSITION_MASK);
				if ((rslts1[i1].flags == kFWD && rslts2[i2].flags == kREV && offset2 >= offset1)
				||  (rslts1[i1].flags == kREV && rslts2[i2].flags == kFWD && offset1 >= offset2))
				{
					unsigned int mm = rslts1[i1].mismatch+rslts2[i2].mismatch;
					if (mm < bestmismatch)
						bestmismatch = mm;
				}
			}
		}
	}
	
	
	// do not gather hits if not worth it.
	if (bestmismatch > kMaxMismatchPerPair)
		goto skipvalidation;
	
	
	for (i1 = 0; i1 < rsltcnt1; i1++)
	{
		unsigned int virus1 = (unsigned int)((rslts1[i1].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT));
		for (i2 = 0; i2 < rsltcnt2; i2++)
		{
			unsigned int virus2 = (unsigned int)((rslts2[i2].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT));
			if ((virus1 == virus2) && ((rslts1[i1].mismatch+rslts2[i2].mismatch) == bestmismatch))
			{
				unsigned int offset1 = (unsigned int)(rslts1[i1].virosauruspos & kPOSITION_MASK);
				unsigned int offset2 = (unsigned int)(rslts2[i2].virosauruspos & kPOSITION_MASK);

				if (rslts1[i1].flags == kFWD && rslts2[i2].flags == kREV)
				{
					if (offset2 >= offset1)
					{
						
						validrslts[validrsltsCnt] = rslts1[i1];
						validrslts[validrsltsCnt].mismatch += rslts2[i2].mismatch; 
						validrslts[validrsltsCnt].taglen = offset2 - offset1 + rslts2[i2].taglen;   
						validrslts[validrsltsCnt].tagid2 = rslts2[i2].tagid;
						validrslts[validrsltsCnt].flags = kFWD;
						validrslts[validrsltsCnt].tag2shift = rslts2[i2].shift;
						validrslts[validrsltsCnt].tag2len = rslts2[i2].tag2len;
					
						if (++validrsltsCnt >= (kMaxHitsPerCPU-8))
						{
							goto skipvalidation;
						}
					}
				}
				else if (rslts1[i1].flags == kREV && rslts2[i2].flags == kFWD)
				{
					if (offset1 >= offset2)
					{
						
						validrslts[validrsltsCnt] = rslts2[i2];
						validrslts[validrsltsCnt].mismatch += rslts1[i1].mismatch; 
						validrslts[validrsltsCnt].taglen = offset1 - offset2 + rslts1[i1].taglen;  
						validrslts[validrsltsCnt].tagid2 = rslts1[i1].tagid;
						validrslts[validrsltsCnt].flags = kREV;
						validrslts[validrsltsCnt].tag2shift = rslts1[i1].shift;
						validrslts[validrsltsCnt].tag2len = rslts1[i1].tag2len;
					
						if (++validrsltsCnt >= (kMaxHitsPerCPU-8))
						{
							goto skipvalidation;
						}
					}
				}
				
				
			}
		}
	}

skipvalidation:;

	
	// keep only smallest insert length of each pair brute force, expecting not so many hits.
	for (i1 = 0; i1 < validrsltsCnt; i1++)
	{
		if (validrslts[i1].mismatch < bestmismatch)
		{
			validrslts[i1].taglen = 0xFFFFFFFF;
			continue;
		}
		if (validrslts[i1].taglen == 0xFFFFFFFF)
			continue;
		unsigned int offset1 = (unsigned int)(validrslts[i1].virosauruspos & kPOSITION_MASK);
		unsigned int virus1 = (unsigned int)(validrslts[i1].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT);
		for (i2 = i1+1; i2 < validrsltsCnt; i2++)
		{
			if (validrslts[i2].taglen == 0xFFFFFFFF)
				continue;
			unsigned int offset2 = (unsigned int)(validrslts[i2].virosauruspos & kPOSITION_MASK);
			unsigned int virus2 = (unsigned int)(validrslts[i2].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT);
			if (offset2 == offset1 && virus2 == virus1) // mask the worst mapping  (e.g. the one with fraglen deviating the most from the most frequent fragment length)
			{
			
				int delta1 = validrslts[i1].taglen - gMostFreqFragLen; if (delta1 < 0) delta1 = -delta1;
				int delta2 = validrslts[i2].taglen - gMostFreqFragLen; if (delta2 < 0) delta2 = -delta2;
				
			
				if (delta2 >= delta1)
					validrslts[i2].taglen = 0xFFFFFFFF;
				else
				{
					validrslts[i1].taglen = 0xFFFFFFFF;
					goto nexti1;
				}
			}
		} 
		nexti1:;
	} 

	int validAltmappingCnt = 0;
	int isMappingOnSingleViruskind = 1; 
	unsigned int currentViruskind = 0xFFFFFFFF;
	unsigned int splitAssignment = 0;
	int multipleVirusHit[kMaxMultGeneHit];

	memset(multipleVirusHit,0,kMaxMultGeneHit*sizeof(int));
	for (i1 = 0; i1 < validrsltsCnt; i1++)
	{
		if (validrslts[i1].taglen == 0xFFFFFFFF)
			continue;
		if (validrslts[i1].mismatch == bestmismatch)
		{
			validAltmappingCnt++;
			unsigned int virus = (unsigned int)(validrslts[i1].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT);
			unsigned int viruskind = ((virus & (kVIRUS_KIND_MASK)) >>  kVIRUS_ENCODING_BITSHIFT);
			if (currentViruskind == 0xFFFFFFFF) 
				currentViruskind = viruskind;
			else
			{
				if (viruskind != currentViruskind)
					isMappingOnSingleViruskind = 0;
			}

			unsigned int i;
			for (i = 0; i< splitAssignment; i++)
			{
				if (virus == multipleVirusHit[i])
					goto alreadyRecorded;
			}
			if (splitAssignment < kMaxMultGeneHit)
			{
				multipleVirusHit[splitAssignment++] = virus;
			}
			else
			{
				printf("WARN: more than %d source genes possible for tag %12lu\t%s\t%s\n",kMaxMultGeneHit,UID,tag1,tag2); fflush(stdout);
				break;
			}
		}
		alreadyRecorded:;
	}

	
	unsigned int reallyValidrsltsCnt = 0;
	for (i1 = 0; i1 < validrsltsCnt; i1++)
	{
		if (validrslts[i1].taglen == 0xFFFFFFFF)
			continue;
		
		validrslts[i1].UID = UID;
		validrslts[i1].splitAssignment = validAltmappingCnt  ;
		validrslts[i1].isMappingOnSingleViruskind = isMappingOnSingleViruskind;
	
		{
			unsigned int k;

			qsort(multipleVirusHit,splitAssignment, sizeof(int), sortMultVirusFunc);
			validrslts[i1].multipleVirusHit = 1;
			int curVirus = multipleVirusHit[0];
			for (k = 0; k < splitAssignment; k++)
			{
				if (multipleVirusHit[k] == 0)
					break;
				if (multipleVirusHit[k] != curVirus)
				{
				
					validrslts[i1].multipleVirusHit++;
				}
				curVirus = multipleVirusHit[k];
			}
		}

		unsigned int virus_chunk = (unsigned int)((validrslts[i1].virosauruspos >> kCHUNK_ENCODING_BITSHIFT) &  (kMAX_VIRUS_CHUNKS-1));

		VIROSAURUSRSLT *mr = gMappingResults[cpu][virus_chunk];
		mr[ gMappingResultsCnt[cpu][virus_chunk].cnt++ ] = validrslts[i1];
		if (gMappingResultsCnt[cpu][virus_chunk].cnt == gMappingResultsCnt[cpu][virus_chunk].maxcnt)
		{

			gMappingResultsCnt[cpu][virus_chunk].maxcnt += kMAPPINGRESULTS_CHUNK;
			gMappingResults[cpu][virus_chunk] = realloc(gMappingResults[cpu][virus_chunk],gMappingResultsCnt[cpu][virus_chunk].maxcnt*sizeof(VIROSAURUSRSLT));
			if (!gMappingResults[cpu][virus_chunk])
			{
				printf("Error: Failed reallocating mapping rslts table %d\n",cpu);
				exit(0);
			}
			else
			{
			}
		}
		reallyValidrsltsCnt++;
	}


	return(reallyValidrsltsCnt);
	
} // ValidatePairs
//---------------------------------------------------------------
int ValidateSingle(unsigned long UID,VIROSAURUSRSLT *rslts, unsigned int rsltcnt,char *tag,char *revtag, char *mate,char *revmate, unsigned int cpu)
{
	unsigned int i1;
	VIROSAURUSRSLT	validrslts[kMaxHitsPerCPU];
	unsigned int validrsltsCnt = 0;
	unsigned int bestmismatch = 999999;

	
	// first pass to determine lowest mismatch
	for (i1 = 0; i1 < rsltcnt; i1++)
	{
		unsigned int mm = rslts[i1].mismatch;
		if (mm < bestmismatch)
			bestmismatch = mm;
	}
	
	
	// do not gather hits if not worth it.
	if (bestmismatch > 8)   // FIXME arbitrary
		goto overflow;
	
	
	for (i1 = 0; i1 < rsltcnt; i1++)
	{
		if (rslts[i1].mismatch == bestmismatch)
		{
			validrslts[validrsltsCnt] = rslts[i1];
			validrslts[validrsltsCnt].taglen = strlen(tag); 
			validrslts[validrsltsCnt].tagid2 = 0xFFFFFFFF; 
			validrslts[validrsltsCnt].flags = rslts[i1].flags; 
			validrslts[validrsltsCnt].tag2shift = 0;
		
			if (++validrsltsCnt >= (kMaxHitsPerCPU-8))
			{
				goto overflow;
			}
		}
	}


overflow:;


	int validAltmappingCnt = 0;
	int isMappingOnSingleViruskind = 1; 
	unsigned int currentViruskind = 0xFFFFFFFF;
	unsigned int splitAssignment = 0;
	int multipleVirusHit[kMaxMultGeneHit];
	memset(multipleVirusHit,0,kMaxMultGeneHit*sizeof(int));
	for (i1 = 0; i1 < validrsltsCnt; i1++)
	{
		validAltmappingCnt++;
		unsigned int virus = (unsigned int)(validrslts[i1].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT);
		unsigned int viruskind = ((virus & (kVIRUS_KIND_MASK)) >>  kVIRUS_ENCODING_BITSHIFT);
		if (currentViruskind == 0xFFFFFFFF) 
			currentViruskind = viruskind;
		else
		{
			if (viruskind != currentViruskind)
				isMappingOnSingleViruskind = 0;
		}

		unsigned int i;
		for (i = 0; i< splitAssignment; i++)
		{
			if (virus == multipleVirusHit[i])
				goto alreadyRecorded;
		}
		if (splitAssignment < kMaxMultGeneHit)
		{
			multipleVirusHit[splitAssignment++] = virus;
		}
		else
		{
			printf("WARN: more than %d source genes possible for tag %12lu\t%s\n",kMaxMultGeneHit,UID,tag); fflush(stdout);
			break;
		}
		alreadyRecorded:;
	}


	
	unsigned int reallyValidrsltsCnt = 0;
	for (i1 = 0; i1 < validrsltsCnt; i1++)
	{
		validrslts[i1].UID = UID;
		validrslts[i1].splitAssignment = validAltmappingCnt  ;
		validrslts[i1].isMappingOnSingleViruskind = isMappingOnSingleViruskind;

		{
			unsigned int k;
	
	
			qsort(multipleVirusHit,splitAssignment, sizeof(int), sortMultVirusFunc);
			validrslts[i1].multipleVirusHit = 1;
			int curVirus = multipleVirusHit[0];
			for (k = 0; k < splitAssignment; k++)
			{
				if (multipleVirusHit[k] == 0)
					break;
				if (multipleVirusHit[k] != curVirus)
				{
				
					validrslts[i1].multipleVirusHit++;
				}
				curVirus = multipleVirusHit[k];
			}
		}
	
		unsigned int virus_chunk = (unsigned int)((validrslts[i1].virosauruspos >> kCHUNK_ENCODING_BITSHIFT) &  (kMAX_VIRUS_CHUNKS-1));

		VIROSAURUSRSLT *mr = gMappingResults[cpu][virus_chunk];
		mr[ gMappingResultsCnt[cpu][virus_chunk].cnt++ ] = validrslts[i1];
		if (gMappingResultsCnt[cpu][virus_chunk].cnt == gMappingResultsCnt[cpu][virus_chunk].maxcnt)
		{

			gMappingResultsCnt[cpu][virus_chunk].maxcnt += kMAPPINGRESULTS_CHUNK;
			gMappingResults[cpu][virus_chunk] = realloc(gMappingResults[cpu][virus_chunk],gMappingResultsCnt[cpu][virus_chunk].maxcnt*sizeof(VIROSAURUSRSLT));
			if (!gMappingResults[cpu][virus_chunk])
			{
				printf("Error: Failed reallocating mapping rslts table %d\n",cpu);
				exit(0);
			}
		}

		reallyValidrsltsCnt++;
	}

	return(reallyValidrsltsCnt);
	
} // ValidateSingle
//---------------------------------------------------------------

void *sendVIROSAURUSfromFastaPairedEnd(void *senddata)
{
	char tag1[kTagSizeStorage];
	char tag2[kTagSizeStorage];
	char revtag1[kTagSizeStorage];
	char revtag2[kTagSizeStorage];
	VIROSAURUSRSLT	tag1rslts[kMaxHitsPerCPU];
	VIROSAURUSRSLT	tag2rslts[kMaxHitsPerCPU];
	unsigned int tag1rsltsCnt;
	unsigned int tag2rsltsCnt;

	unsigned int cpu = ((THREAD_DATA*)senddata)->cpu;
	unsigned char *queryTag512bits = NULL;
	unsigned int bufpos = cpu * kINPUTBUFFER_MAXCNT_PER_CPU;
	unsigned long UID;
	int err;
	gCrappy[cpu] = 0;
	gUnmappedCnt[cpu] = 0;

	err  = posix_memalign((void **)&queryTag512bits, getpagesize(),	(unsigned long)(64)); 
	if (err)
		goto bail;

	// ---------  infinite loop  --------------------------------
	do
	{

		// find non-empty slot
		while(inputbuffer[bufpos*kMaxInputLineLen+kMaxInputLineLen-1] == 0)
		{
			if (inputbuffer[cpu*kINPUTBUFFER_MAXCNT_PER_CPU+kMaxInputLineLen-2] == 99)
				goto bail;
			
			bufpos++;
			if (bufpos == (cpu+1) * kINPUTBUFFER_MAXCNT_PER_CPU)
				bufpos = cpu * kINPUTBUFFER_MAXCNT_PER_CPU;
			
		};
		sscanf(&inputbuffer[bufpos*kMaxInputLineLen],"%lu\t%s\t%s\n",&UID,&tag1[0],&tag2[0]);
		inputbuffer[bufpos*kMaxInputLineLen+kMaxInputLineLen-1]--;
		
		// CHECK LENGTH
		{
			int  k;
			for (k=0; k<gMinTagLen;k++)
			{
				if (tag1[k] == 0 || (gSingleEnd==0 && tag2[k] == 0))
				{
					gCrappy[cpu]++;
					goto skipCrappy;
				}
			}
		}

		
		// filter crap
		{

			char *hit;
			unsigned int hitCnt;
			#define kRepeatFilterThreshold 9
			
			hit = &tag1[0]; hitCnt = 0;
			while (hit = strstr(hit,"TTAGGG")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

			hit = &tag1[0]; hitCnt = 0;
			while (hit = strstr(hit,"CCCTAA")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

			hit = &tag1[0]; hitCnt = 0;
			while (hit = strstr(hit,"CTTCTC")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

			hit = &tag1[0]; hitCnt = 0;
			while (hit = strstr(hit,"CCTTCT")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

			hit = &tag1[0]; hitCnt = 0;
			while (hit = strstr(hit,"GAGAGA")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

			hit = &tag1[0]; hitCnt = 0;
			while (hit = strstr(hit,"CACACA")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

			hit = &tag1[0]; hitCnt = 0;
			while (hit = strstr(hit,"TTTTTT")) { hitCnt++; hit+=6; } ; if (hitCnt >= 10) { gCrappy[cpu]++; goto skipCrappy; }

			hit = &tag1[0]; hitCnt = 0;
			while (hit = strstr(hit,"AAAAAA")) { hitCnt++; hit+=6; } ; if (hitCnt >= 10) { gCrappy[cpu]++; goto skipCrappy; }

			if (gSingleEnd == 0)
			{
				hit = &tag2[0]; hitCnt = 0;
				while (hit = strstr(hit,"TTAGGG")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

				hit = &tag2[0]; hitCnt = 0;
				while (hit = strstr(hit,"CCCTAA")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

				hit = &tag2[0]; hitCnt = 0;
				while (hit = strstr(hit,"CTTCTC")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

				hit = &tag2[0]; hitCnt = 0;
				while (hit = strstr(hit,"CCTTCT")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

				hit = &tag2[0]; hitCnt = 0;
				while (hit = strstr(hit,"GAGAGA")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

				hit = &tag2[0]; hitCnt = 0;
				while (hit = strstr(hit,"CACACA")) { hitCnt++; hit+=6; } ; if (hitCnt >= kRepeatFilterThreshold) { gCrappy[cpu]++; goto skipCrappy; }

				hit = &tag2[0]; hitCnt = 0;
				while (hit = strstr(hit,"TTTTTT")) { hitCnt++; hit+=6; } ; if (hitCnt >= 10) { gCrappy[cpu]++; goto skipCrappy; }

				hit = &tag2[0]; hitCnt = 0;
				while (hit = strstr(hit,"AAAAAA")) { hitCnt++; hit+=6; } ; if (hitCnt >= 10) { gCrappy[cpu]++; goto skipCrappy; }
			}
		}
		
		unsigned int minMismatch1;
		tag1rsltsCnt = DoSearchVirosaurus2(&tag1[0],&revtag1[0],&queryTag512bits[0],&tag1rslts[0],&minMismatch1);


		tag1rsltsCnt = VirosaurusCheckHits2(&tag1rslts[0] ,tag1rsltsCnt ,indexLeft ,&tag1[0],&revtag1[0],minMismatch1);


		if (gSingleEnd)
			tag2rsltsCnt = 0;
		else
		{
	
			unsigned int minMismatch2;
			tag2rsltsCnt = DoSearchVirosaurus2(&tag2[0],&revtag2[0],&queryTag512bits[0],&tag2rslts[0],&minMismatch2);
			tag2rsltsCnt = VirosaurusCheckHits2(&tag2rslts[0] ,tag2rsltsCnt ,indexLeft ,&tag2[0],&revtag2[0],minMismatch2);
		}
		if ((tag1rsltsCnt == 0) && (tag2rsltsCnt > 0))
		{
		
			if (ValidateSingle(UID,&tag2rslts[0],tag2rsltsCnt,tag2,&revtag2[0],tag1,&revtag1[0],cpu) == 0)
				goto uglyunmapped;
			else
				goto skipCrappy;
		}
		else if ((tag2rsltsCnt == 0) && (tag1rsltsCnt > 0))
		{
		
			if (ValidateSingle(UID,&tag1rslts[0],tag1rsltsCnt,tag1,&revtag1[0],tag2,&revtag2[0],cpu) == 0)
				goto uglyunmapped;
			else
				goto skipCrappy;
		}


		if (tag1rsltsCnt > 0 && tag2rsltsCnt > 0)
		{
			if (ValidatePairs(rsltfile[cpu],UID,&tag1rslts[0] ,tag1rsltsCnt,&tag2rslts[0] ,tag2rsltsCnt,tag1,tag2,&revtag1[0],&revtag2[0],cpu) == 0)
			{
				
				if (ValidateSingle(UID,&tag2rslts[0],tag2rsltsCnt,tag2,&revtag2[0],tag1,&revtag1[0],cpu) != 0)
					goto skipCrappy;
				if (ValidateSingle(UID,&tag1rslts[0],tag1rsltsCnt,tag1,&revtag1[0],tag2,&revtag2[0],cpu) == 0)
					goto skipCrappy;
				goto uglyunmapped;
			}
		}
		else
		{
uglyunmapped:;
			
			{
				gUnmappedCnt[cpu]++;
			}
		}
	skipCrappy:;
		gProcessedCnt[cpu]++;
	
		bufpos++;
		if (bufpos == (cpu+1) * kINPUTBUFFER_MAXCNT_PER_CPU)
			bufpos = cpu * kINPUTBUFFER_MAXCNT_PER_CPU;

	} while(1);


	

bail:;

	
	if (queryTag512bits)
		free(queryTag512bits);

	
	return(0);
	
} // sendVIROSAURUSfromFastaPairedEnd
//---------------------------------------------------------------

void *readingSingleEnd(void *fn)
{

	unsigned int bufpos = 0;
	unsigned long volatile sent = 0;
	gMaxUID = 0;
	char R1[512];

	if ( ((char*)fn)[0] == 0)
		return(0);
		
		
	printf("Opening jobfile %s\n",(char*)fn); fflush(stdout);
	FILE *F = fopen((char*)fn,"r");
	if (!F)
	{
		printf("FATAL Cannot open input stream %s\n",(char*)fn); fflush(stdout);
		exit(1);
	}
	
	// ---------  infinite loop  --------------------------------
	do
	{
		
		// wait for empty slot
		while(inputbuffer[bufpos*kMaxInputLineLen+kMaxInputLineLen-1] != 0)
		{
			bufpos++;
			if (bufpos == gThreadCnt*kINPUTBUFFER_MAXCNT_PER_CPU)
				bufpos = 0;
		}


		// skip read header
		fgets(&R1[0],511,F);
		fgets(&R1[0],511,F);
		
		if (gFilterAdapter[0])
		{
			
			char *stripPtr1 = strstr(R1,gFilterAdapter);
			if (stripPtr1)
			{
				// avoid zero length string.
				if (stripPtr1 == &R1[0])
					stripPtr1++;
				*stripPtr1 = '\0';
			}
		}
		sprintf(&inputbuffer[bufpos*kMaxInputLineLen],"%lu\t%s\n",++gMaxUID,R1);

		// skip + and quality score
		fgets(&R1[0],511,F);
		fgets(&R1[0],511,F);
		if (feof(F))
			goto readingEnd;
		
		
		inputbuffer[bufpos*kMaxInputLineLen+kMaxInputLineLen-1]++;

		sent++;
		bufpos++;
		if (bufpos == gThreadCnt*kINPUTBUFFER_MAXCNT_PER_CPU)
			bufpos = 0;

		if ((sent % 20000) == 0)
		{
			gettimeofday(&te, NULL);
			int seconds = ((int)te.tv_sec-(int)ts.tv_sec);
			fprintf(stderr,"\033[s\033[34;1fSent cnt=%10lu %8.1f/s\033[u",sent,sent/((float)seconds));
			fflush(stderr);
		}

	} while(1);
	// ---- end of endless loop.


readingEnd:;

	fclose(F);
	printf("\nDone with file %s ; sent %lu pairs; signalling end of jobs\n",(char*)fn,sent);
	fflush(stdout);
	sleep(1);

	
	// wait until all processed.
	unsigned long volatile totProcessed;
	unsigned long crappy = 0;
	unsigned long unmapped = 0;
	while(1)
	{
		totProcessed = 0;
		for (bufpos = 0; bufpos < gThreadCnt; bufpos++)
		{
			totProcessed += gProcessedCnt[bufpos];
		}
		if (totProcessed == sent)
			break;

		fprintf(stderr,"Waiting (processed cnt=%10lu / %10lu )\n",totProcessed,sent); fflush(stderr);
		sleep(1);
	};
	
	for (bufpos = 0; bufpos < gThreadCnt; bufpos++)
	{
		inputbuffer[bufpos*kINPUTBUFFER_MAXCNT_PER_CPU+kMaxInputLineLen-2]=99;
		crappy += gCrappy[bufpos];
		unmapped += gUnmappedCnt[bufpos];
	}

	printf("Filtered %lu crappy pairs out of %lu pairs. Unfiltered but unmapped: %lu\n",crappy,sent,unmapped);
	fflush(stdout);


	return(0);

} // readingSingleEnd
//---------------------------------------------------------------

void *readingPairedEnd(void *fn)
{

	unsigned int bufpos = 0;
	unsigned long volatile sent = 0;
	gMaxUID = 0;
	char R1[512];
	char R2[512];

	if ( ((char*)fn)[0] == 0)
		return(0);
		
		
	printf("Opening jobfile %s\n",(char*)fn); fflush(stdout);
	FILE *F = fopen((char*)fn,"r");
	if (!F)
	{
		printf("FATAL Cannot open input stream %s\n",(char*)fn); fflush(stdout);
		exit(1);
	}
	FILE *F2 = fopen(gInputFile2,"r");
	if (!F2)
	{
		printf("FATAL Cannot open input stream %s\n",gInputFile2); fflush(stdout);
		fclose(F);
		exit(1);
	}
	
	do
	{
		// wait for empty slot
		while(inputbuffer[bufpos*kMaxInputLineLen+kMaxInputLineLen-1] != 0)
		{
			bufpos++;
			if (bufpos == gThreadCnt*kINPUTBUFFER_MAXCNT_PER_CPU)
				bufpos = 0;
		}
		
		// skip read header
		fgets(&R1[0],511,F);
		fgets(&R2[0],511,F2);

		fgets(&R1[0],511,F);
		fgets(&R2[0],511,F2);
		
		if (gFilterAdapter[0])
		{
			
			char *stripPtr1 = strstr(R1,gFilterAdapter);
			if (stripPtr1)
			{
				
				if (stripPtr1 == &R1[0])
					stripPtr1++;
				*stripPtr1 = '\0';
			}
			char *stripPtr2 = strstr(R2,gFilterAdapter);
			if (stripPtr2)
			{
				
				if (stripPtr2 == &R2[0])
					stripPtr2++;
				*stripPtr2 = '\0';
			}
		}
		sprintf(&inputbuffer[bufpos*kMaxInputLineLen],"%lu\t%s\t%s\n",++gMaxUID,R1,R2);

		
		// skip + and quality score
		fgets(&R1[0],511,F);
		fgets(&R1[0],511,F);
		fgets(&R2[0],511,F2);
		fgets(&R2[0],511,F2);
		if (feof(F))
			goto readingEnd;
		
		
		inputbuffer[bufpos*kMaxInputLineLen+kMaxInputLineLen-1]++;

		sent++;
		bufpos++;
		if (bufpos == gThreadCnt*kINPUTBUFFER_MAXCNT_PER_CPU)
			bufpos = 0;

		if ((sent % 20000) == 0)
		{
			gettimeofday(&te, NULL);
			int seconds = ((int)te.tv_sec-(int)ts.tv_sec);
			fprintf(stderr,"\033[s\033[34;1fSent cnt=%10lu %8.1f/s\033[u",sent,sent/((float)seconds));
			fflush(stderr);
		}


	} while(1);
	// ---- end of endless loop.

	

readingEnd:;

	fclose(F2);
	fclose(F);
	printf("\nDone with file %s ; sent %lu pairs; signalling end of jobs\n",(char*)fn,sent);
	fflush(stdout);
	sleep(1);

	
	// wait until all processed.
	unsigned long volatile totProcessed;
	unsigned long crappy = 0;
	unsigned long unmapped = 0;
	while(1)
	{
		totProcessed = 0;
		for (bufpos = 0; bufpos < gThreadCnt; bufpos++)
		{
			totProcessed += gProcessedCnt[bufpos];
		}
		if (totProcessed == sent)
			break;

		fprintf(stderr,"Waiting (processed cnt=%10lu / %10lu )\n",totProcessed,sent); fflush(stderr);
		sleep(1);
	};
	
	for (bufpos = 0; bufpos < gThreadCnt; bufpos++)
	{
		inputbuffer[bufpos*kINPUTBUFFER_MAXCNT_PER_CPU+kMaxInputLineLen-2]=99;
		crappy += gCrappy[bufpos];
		unmapped += gUnmappedCnt[bufpos];
	}

	printf("Filtered %lu crappy pairs out of %lu pairs. Unfiltered but unmapped: %lu\n",crappy,sent,unmapped);
	fflush(stdout);

	return(0);

	
} // readingPairedEnd
//---------------------------------------------------------------

typedef struct pos_index
{
    unsigned long pos;
    unsigned int uid;
    unsigned int idx;
	unsigned int srccpu;
} t_pos_index;
 
static int compare_data(const void *a, const void *b);

static int compare_data (const void *a, const void *b)
{
     
	t_pos_index *struct_a = (t_pos_index *) a;
	t_pos_index *struct_b = (t_pos_index *) b;

	if(struct_a->pos < struct_b->pos)
		return -1;
	if(struct_a->pos == struct_b->pos)
	{

		if(struct_a->uid < struct_b->uid)
			return -1;
		else
			return 1;
	}
	return 1;
      
}  
//---------------------------------------------------------------

void *sortMappingResults(void *senddata)
{

	unsigned int chunk = ((THREAD_DATA*)senddata)->cpu;
	unsigned int filterVirusID = ((THREAD_DATA*)senddata)-> filterVirusID;
	unsigned int cpu;
	unsigned int i;
	unsigned int datacnt;
	t_pos_index *dataA = NULL;
	FILE *f = NULL;
	char rsltfilename[kMaxFileName+32];

	unsigned int chunkTot = 0;
	for (cpu = 0; cpu < gThreadCnt; cpu++)
		chunkTot += gMappingResultsCnt[cpu][chunk].cnt;

	if (chunkTot == 0)
		goto empty;
	
	dataA = malloc(chunkTot*sizeof(t_pos_index));
	if (!dataA)
	{
		printf("out of mem\n");
		exit(1);
	}
	unsigned char *multipleMappingHitFlags = calloc((gMaxUID+1),sizeof(char));
	if (!multipleMappingHitFlags)
	{
		printf("out of mem for flags\n");
		exit(1);
	}

	datacnt = 0;
	for (cpu = 0; cpu < gThreadCnt; cpu++)
	{
		VIROSAURUSRSLT *mr = gMappingResults[ cpu ][ chunk ];
		for (i = 0; i < gMappingResultsCnt[cpu][chunk].cnt; i++)
		{
			unsigned int virus = (unsigned int)((mr[i].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT));
			unsigned int offset = (unsigned int)(mr[i].virosauruspos & kPOSITION_MASK)  -mr[i].shift  ;


			if (mr[i].tagid2 == 0xFFFFFFFF) 
			{
			}
			else
			{
				unsigned long start2;
				start2 = &virosaurusDB[virosaurusDBidx[mr[i].tagid2]-mr[i].tag2shift]  - &virosaurusDB[virosaurusDBidx[mr[i].tagid]-mr[i].shift];
				mr[i].taglen = (unsigned int)start2 + mr[i].tag2len;
			}


			dataA[datacnt].pos = (unsigned long)(virus);  
			dataA[datacnt].pos <<= 22;
			dataA[datacnt].pos += (unsigned long)offset & 0x3FFFFF;
			
			dataA[datacnt].pos <<= 12;
			dataA[datacnt].pos += (mr[i].taglen & 0x0FFF);
			dataA[datacnt].pos <<= 4; 
			dataA[datacnt].pos += mr[i].mismatch;  
			dataA[datacnt].uid = mr[i].UID;
			dataA[datacnt].idx = i;
			dataA[datacnt].srccpu = cpu;
			datacnt++;

		}
	}


	qsort(dataA, chunkTot, sizeof(t_pos_index), compare_data);


	sprintf(rsltfilename,"%s/chunk%d.rslts",gOutputDir,1+chunk);
	f = fopen(rsltfilename,"w");
	if (f == NULL)
	{
		printf("cannot write file %s\n",rsltfilename);
		exit(1);
	}
	
	
	VIROSAURUSRSLT *lastwritten_ptr = NULL;
	unsigned int currentVirusFull = 0xFFFFFFFF;
	unsigned int currentViruskind = 0xFFFFFFFF;
	unsigned int nrcnt = 0;
	char multMapFlag;	
	

	unsigned int virussubtype;
	unsigned int offset;
	unsigned int viruskind;
	
		
		for (i = 0; i < chunkTot; i++)
		{

			VIROSAURUSRSLT *mr = gMappingResults[ dataA[i].srccpu ][ chunk ];
			unsigned int i1 = dataA[i].idx;
			unsigned int fraglen;

			unsigned int virusfull = (unsigned int)((mr[i1].virosauruspos >> kVIRUS_TYPEandSUBTYPE_BITSHIFT));
			virussubtype = virusfull & (kMAX_VIRUS_SUBTYPES-1);
			viruskind = ((virusfull & (kVIRUS_KIND_MASK)) >>  kVIRUS_ENCODING_BITSHIFT);
			offset = (unsigned int)(mr[i1].virosauruspos & kPOSITION_MASK)   -mr[i1].shift  ;
			fraglen = mr[i1].taglen;

				if ((lastwritten_ptr == NULL) ||  mr[i1].virosauruspos != lastwritten_ptr->virosauruspos ||  mr[i1].taglen != lastwritten_ptr->taglen)
				{
				
				
					if (virusfull == currentVirusFull)
					{
						if ((multipleMappingHitFlags[mr[i1].UID] & 0x01) != 0)
						{
								multMapFlag = 'i';
								goto pairedEndTagAlreadyWritten;
						}
						multipleMappingHitFlags[mr[i1].UID] |= 0x01;
					}
					else
					{
						
						if (viruskind == currentViruskind)
						{
							
							unsigned int k;
							for (k=0; k< (gMaxUID+1); k++)
								multipleMappingHitFlags[k] &= 0xFE;
						}
						else
						{
							memset(multipleMappingHitFlags,0,(gMaxUID+1)*sizeof(char));
							currentViruskind = viruskind;
						}
						currentVirusFull = virusfull;
					}
					nrcnt++;
					multMapFlag = 'k'; 

		pairedEndTagAlreadyWritten:;
					lastwritten_ptr = &mr[i1];

					if (multMapFlag == 'k' && viruskind < filterVirusID)
						fprintf(f,"%u,%u\t%d\t%u\t%u\t%u\t%u\t%d\n",viruskind,virussubtype,offset,fraglen,mr[i1].mismatch,mr[i1].splitAssignment,mr[i1].multipleVirusHit,mr[i1].isMappingOnSingleViruskind);


					{

						if (virussubtype < kMAX_VIRUS_SUBTYPES)
						{
							gVirusCnt[viruskind*kMAX_VIRUS_SUBTYPES+virussubtype].virus += 1.0/mr[i1].splitAssignment ;
							if (mr[i1].isMappingOnSingleViruskind && mr[i1].mismatch <= gMaxMismatchForUnambiguous && mr[i1].taglen >= gMinFragLen)
								gVirusCnt[viruskind*kMAX_VIRUS_SUBTYPES+virussubtype].virusUniqueGood += 1.0/mr[i1].splitAssignment;
							gVirusCnt[viruskind*kMAX_VIRUS_SUBTYPES+virussubtype].nameidx = mr[i1].tagid; 
						}
						else
						{
							printf("ERR virus cluster overflow: %u\n",virussubtype);
							fflush(stdout);
						}
					}
					
					multipleMappingHitFlags[mr[i1].UID] |= 0x02; 
					multipleMappingHitFlags[mr[i1].UID] |= 0x01; 
				}
				else
				{
				}
		}

	if (nrcnt > 0)
	{	printf("wrote %12u non-redundant mappings out of %12u mappings to chunk%3u [%7.3f %%] file %s\n",nrcnt,chunkTot,1+chunk,100.0*(float)nrcnt/chunkTot,rsltfilename); fflush(stdout); }

	fclose(f);
	free(dataA);
	free(multipleMappingHitFlags);

empty:;

	for (cpu = 0; cpu < gThreadCnt; cpu++)
		free(gMappingResults[cpu][chunk]);

	return(0);

} // sortMappingResults

//---------------------------------------------------------------
int main (int argc, char **argv)
{
	HEADER hdr1;

	int c;
	int err = 0;
	char virosaurusfile[kMaxFileName];
	char rsltfilename[kMaxFileName+32];
	size_t memSize;
	FILE *f;
	unsigned int cpu;
	unsigned int chr;
	pthread_t readingThread;
	pthread_t sortingthread[kMAX_VIRUS_CHUNKS];
	pthread_t thread[kMAX_TAG_SNDRCV_THREAD_cnt];
	THREAD_DATA  senddata[kMAX_TAG_SNDRCV_THREAD_cnt];

	size_t ps = (size_t)getpagesize();
	virosaurusfile[0] = 0;
	gInputFile[0] = 0;
	gInputFile2[0] = 0;
	gFilterAdapter[0] = 0;
	strcpy(gOutputDir,"./");
	gMinTagLen = 36;
	gMinFragLen = 101;
	gMaxMismatchForUnambiguous = 4;
	gThreadCnt = 1;

	
	opterr = 0;
	while ((c = getopt (argc, argv, "1:2:d:t:f:L:F:M:")) != -1)
	switch (c)
	{      

	  case '1':
			strcpy(gInputFile,optarg);
        break;
	  case '2':
			strcpy(gInputFile2,optarg);
        break;

	  case 'd':
			strcpy(gOutputDir,optarg);
			break;

	  case 't':
			sscanf(optarg,"%u",&gThreadCnt);
			if (gThreadCnt > kMAX_TAG_SNDRCV_THREAD_cnt)
				gThreadCnt = kMAX_TAG_SNDRCV_THREAD_cnt;
        break;

	  case 'f':
			strcpy(gFilterAdapter,optarg);
        break;
		
		case 'L':
			sscanf(optarg,"%u",&gMinTagLen);
		break;

		case 'F':
			sscanf(optarg,"%u",&gMinFragLen);
		break;

		case 'M':
			sscanf(optarg,"%u",&gMaxMismatchForUnambiguous);
		break;
		
	}

	if ((gInputFile[0] == 0) || gInputFile2[0] == 0)
	{
		fprintf(stderr,"Usage:\n");

		fprintf(stderr,"virus_match   -1 read1 -2 read2 [-f filterSequence][-L minTagLen][-F minFragLen][-M maxMismatch][-t threadCnt][-d output_directory]\n\n");
	
		fprintf(stderr,"  -1             paired-end fastq file of read1\n");
		fprintf(stderr,"  -2             paired-end fastq file of read2\n");

		fprintf(stderr,"  -f sequence    filter adapter sequence\n");
		fprintf(stderr,"                 DNA: illumina nextera   = CTGTCTCTTATACACA TCT\n");
		fprintf(stderr,"                 RNA: truseq riboprofile = AGATCGGAAGAGC ACACGTCT\n");
	
		fprintf(stderr,"  -L minTagLen   discard any pair if tag1 or tag2 is below this value [default: %u nt]\n",gMinTagLen);
		fprintf(stderr,"  -F minFragLen  minimal fragment length to record a hit as unambiguous in the stats [default: %u nt]\n",gMinFragLen);
		fprintf(stderr,"  -M maxMismatch maximal mismatch count for fragment to record a hit as unambiguous in the stats [default: %u nt]\n",gMaxMismatchForUnambiguous);
		fprintf(stderr,"  -d             output directory (must exist default: ./)\n");
		fprintf(stderr,"  -t threadCnt   number of concurrent threads (max=%d)\n",kMAX_TAG_SNDRCV_THREAD_cnt);
	

		fprintf(stderr,"\nExample:\nOUTDIR=/tmp/virusTEST; mkdir -p $OUTDIR; virus_match_AVX512 -d $OUTDIR -f AGATCGGAAGAGC -L 100 -M 12  -1 <(gzip -dc test_R1.fastq.gz) -2 <(gzip -dc test_R2.fastq.gz) -t 64 > $OUTDIR/log\n\n");

		return(0);
	}


	char szTmp[64];
	sprintf(szTmp, "/proc/%d/exe", getpid());
	int bytes = readlink(szTmp, virosaurusfile, kMaxFileName);
	strcpy(&virosaurusfile[bytes-18],"20180906.bin");


	if (gInputFile2[0] == 0)
		gSingleEnd = 1;


	{	printf("VIRUS_MATCH [HUG VERSION=1.0]: N.Guex & C.Iseli\nVIRAL_GENOME=%s: P.Lemercier & A.Gleize\nFilterAdapter=%s; MinTagLen=%u; MinFragLen=%u; MaxMismatch=%u\n",virosaurusfile,gFilterAdapter,gMinTagLen,gMinFragLen,gMaxMismatchForUnambiguous); fflush(stdout);}

	
	// ---------  start loading input buffer --------------------------------

	err = posix_memalign((void **)&inputbuffer,  getpagesize(), kINPUTBUFFER_DATA_SIZE);
	if (err != 0)
	{
		printf("failed allocating memory for input buffer.\n");
		goto bail;
	}
	memset(inputbuffer,0,kINPUTBUFFER_DATA_SIZE);
	if (gSingleEnd)
	{
		if (pthread_create (&readingThread, NULL, &readingSingleEnd, (void*)&gInputFile))
		{
			printf("Error: Failed creating readingThread thread\n");
			goto bail;
		}
	}
	else
	{
		if (pthread_create (&readingThread, NULL, &readingPairedEnd, (void*)&gInputFile))
		{
			printf("Error: Failed creating readingThread thread\n");
			goto bail;
		}
	}
	

	if (openPackedTags(virosaurusfile, &hdr1) != 0)
		goto bailout;

	
	gettimeofday(&ts, NULL);

	memSize = ((unsigned long)hdr1.tagcnt*k256bits + ps) & (~(ps-1));
	

	err += posix_memalign((void **)&sortedtagsLeft, ps,	memSize);

	memSize = (hdr1.hashSize*sizeof(HASH) + ps) & (~(ps-1));
	err += posix_memalign((void **)&hashLeft, ps,	memSize);

	if (err)
		goto bailout;

	indexLeft = malloc((hdr1.tagcnt+15)*sizeof(unsigned int));
	if (!indexLeft)
		goto bailout;

	tagNames1 = malloc((long)hdr1.tagcnt*(long)kTagNameMaxLen*sizeof(char));
	if (!tagNames1)
		goto bailout;

	tagUID1 = malloc((long)hdr1.tagcnt*sizeof(unsigned long));
	if (!tagUID1)
		goto bailout;

	virosaurusDBidx = malloc((long)hdr1.tagcnt*sizeof(unsigned long));
	if (!virosaurusDBidx)
		goto bailout;

	virosaurusDB = malloc((long)hdr1.virosaurusDBsize*sizeof(char));
	if (!virosaurusDB)
		goto bailout;


	virusshortnames = malloc((hdr1.maxVirusID+1) * kShortNameLen*sizeof(char));
	if (!virusshortnames)
		goto bailout;


	if (hdr1.filterVirusID == 0)
		hdr1.filterVirusID = hdr1.maxVirusID + 1;
	printf("LOADING DB (version %u) maxVirusID=%u filterFromVirusID=%u\n",hdr1.version,hdr1.maxVirusID,hdr1.filterVirusID); fflush(stdout);

	if (loadVirosaurus(tagUID1,tagNames1,sortedtagsLeft,indexLeft,hashLeft,hdr1.tagcnt,virosaurusDB,virosaurusDBidx,virosaurusfile) != 0)
		goto bailout;


	int r;
	for (r = hdr1.tagcnt; r < hdr1.tagcnt+15; r++)
	{
		indexLeft[ r ] = 0;
	}

	gettimeofday(&te, NULL);

	printf("%8d seconds:Loaded virus database\n",((int)te.tv_sec-(int)ts.tv_sec)); fflush(stdout);

		printf("INIT\n"); fflush(stdout);

	initRevCompTable((unsigned char *)gIUPACrevcomp);


	err = posix_memalign((void **)&mem512,    64, 16*sizeof(unsigned int));
	err  += posix_memalign((void **)&masks512, getpagesize(),	kTagSizeStorage*(unsigned long)(64));
	if (err)
		goto bail;

	initPackNT1bitTable(gNT1bit);
	initPackRevNT1bitTable(grevNT1bit);
	init512bitsVectors();
	init512bitsMasks();
	

	gettimeofday(&ts, NULL);

	
	/* ----------------------- launch tag sender threads ----------------------*/

	gVirusCnt = calloc((hdr1.maxVirusID+1)*kMAX_VIRUS_SUBTYPES,sizeof(VIRUSCNT));
	if (gVirusCnt == NULL)
	{
		printf("not enough RAM\n");
		goto bail;
	}

	
	memhash = hashLeft;
	packedtags = sortedtagsLeft;
	
	for (cpu = 0; cpu < gThreadCnt; cpu++) 
	{	
		
		for (chr = 0; chr < kMAX_VIRUS_CHUNKS; chr++)
		{
			gMappingResultsCnt[cpu][chr].cnt = 0;
	
				gMappingResultsCnt[cpu][chr].maxcnt = kMAPPINGRESULTS_CHUNK;
			gMappingResults[cpu][chr] = malloc(gMappingResultsCnt[cpu][chr].maxcnt*sizeof(VIROSAURUSRSLT));
			if (!gMappingResults[cpu][chr])
			{
				printf("Error: Failed creating mapping rslts table %d\n",cpu);
				goto bail;
			}
		}
		senddata[cpu].cpu = cpu;
		gProcessedCnt[cpu] = 0;
		{
			if (pthread_create (&thread[cpu], NULL, &sendVIROSAURUSfromFastaPairedEnd, (void*)&senddata[cpu]))
			{
				printf("Error: Failed creating sendTags thread %d\n",cpu);
				goto bail;
			}
		}
		
	}


	sprintf(rsltfilename,"%s/chunk0.rslts",gOutputDir);
	f = fopen(rsltfilename,"w");
	if (f == NULL)
	{
		printf("cannot write file %s\n",rsltfilename);
		exit(1);
	}
	fprintf(f,"virus\tposition\tfraglen\tmismatch\tbestMapCnt\tdistinctVirusesSubtypesCnt\tsingleVirusFlag\n");
	fclose(f);


	for (cpu = 0; cpu < gThreadCnt; cpu++)
	{
		if (pthread_join (thread[cpu], NULL))
		{
			printf("Error: Failed pthread_join %d\n",cpu);
			goto bail;
		}
	}

	gettimeofday(&te, NULL);
	printf("%8d seconds:Collecting results\n",((int)te.tv_sec-(int)ts.tv_sec)); fflush(stdout);
	for (chr = 0; chr < kMAX_VIRUS_CHUNKS; chr++)
	{
		senddata[chr].cpu = chr;
		senddata[chr].filterVirusID = hdr1.filterVirusID;
		if (pthread_create (&sortingthread[chr], NULL, &sortMappingResults, (void*)&senddata[chr]))
		{
			printf("Error: Failed creating sorting thread %d\n",chr);
			goto bail;
		}
	}


	for (chr = 0; chr < kMAX_VIRUS_CHUNKS; chr++)
	{
		if (pthread_join (sortingthread[chr], NULL))
		{
			printf("Error: Failed pthread_join %d\n",chr);
			goto bail;
		}
	}


	// ---------  END  --------------------------------

	
bail:;

	if (pthread_join (readingThread, NULL))
		printf("Error: Failed pthread_join of readingThread\n");
	

	if (gVirusCnt)
	{
		FILE *F1 = NULL;
		FILE *F2 = NULL;
		char fn[kMaxFileName+32];

		gettimeofday(&te, NULL);
		printf("%8d seconds:Creating Summary\n",((int)te.tv_sec-(int)ts.tv_sec)); fflush(stdout);

		sprintf(fn,"%s/summary_virus.tsv",gOutputDir);
		F1 = fopen(fn,"w");
		if (F1 == NULL)
		{
			printf("cannot write file %s\n",fn);
			goto skip;
		}
		sprintf(fn,"%s/summary_virus-cluster.tsv",gOutputDir);
		F2 = fopen(fn,"w");
		if (F2 == NULL)
		{
			printf("cannot write file %s\n",fn);
			goto skip;
		}

		fprintf(F1,"virus\tID\tsignal\tunambiguous_signal\n");
		fprintf(F2,"virus\tcluster\tID\tsignal\tunambiguous_signal\n");
		for (chr = 0; chr <= hdr1.maxVirusID; chr++)
		{
			
			unsigned int k;
			float sum = 0.0;
			float sumUniqueGood = 0.0;
			for (k = 1; k < kMAX_VIRUS_SUBTYPES; k++)
			{
				
				if (gVirusCnt[chr*kMAX_VIRUS_SUBTYPES+k].virus > 0.0)
				{
					fprintf(F2,"%u\t%u\t%s\t%f\t%f\n",chr,k,&tagNames1[gVirusCnt[chr*kMAX_VIRUS_SUBTYPES+k].nameidx*kTagNameMaxLen],gVirusCnt[chr*kMAX_VIRUS_SUBTYPES+k].virus,gVirusCnt[chr*kMAX_VIRUS_SUBTYPES+k].virusUniqueGood);
				
					sum += gVirusCnt[chr*kMAX_VIRUS_SUBTYPES+k].virus;
					sumUniqueGood += gVirusCnt[chr*kMAX_VIRUS_SUBTYPES+k].virusUniqueGood;
				}
			}
			if (sum > 0.0)
			{
				fprintf(F1,"%u\t%s\t%f\t%f\n",chr,&virusshortnames[chr*kShortNameLen],sum,sumUniqueGood);
			}
		}

skip:;

		if (F1)
			fclose(F1);
		if (F2)
			fclose(F2);

		free(gVirusCnt);
	}
bailout:;
	
	if (tagNames1)
		free(tagNames1);

	if (virosaurusDBidx)
		free(virosaurusDBidx);

	if (virosaurusDB)
		free(virosaurusDB);

	if (tagUID1)
		free(tagUID1);
	
	if (indexLeft)
		free(indexLeft);

	if (hashLeft)
		free(hashLeft);

	if (inputbuffer)
		free(inputbuffer);

	if (mem512)
		free(mem512);

	if (masks512)
		free(masks512);

	if (virusshortnames)
		free(virusshortnames);

	printf("DONE.\n"); fflush(stdout);

	return(err);
	
} // main
//---------------------------------------------------------------
// vim: tabstop=2 shiftwidth=2
