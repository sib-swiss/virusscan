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


#define k512bits 64
#define k256bits 32
#define k128bits 16


#define kHashEncodingPos 0

#define kHashSize 0x00008000
#define kHashMask 0x00007FFF
#define kHashShift  0

#define kTagNameMaxLen 64L
#define kShortNameLen 64
#define kDatabaseTagStorageSize 128  /* !!!! Warning !!!!  assembly code hardwired for this size */

#define kPOSITION_MASK 0x00000000FFFFFFFF
#define kMAX_DISTINCT_VIRUSES 0x2000
#define kMAX_VIRUS_SUBTYPES 0x2000
#define kVIRUS_KIND_MASK 0x3FFFFFF
#define kVIRUS_ENCODING_BITSHIFT 13
#define kCHUNK_ENCODING_BITSHIFT 58
#define kVIRUS_TYPEandSUBTYPE_BITSHIFT 32
#define kMAX_VIRUS_CHUNKS 256

/*

xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
[ chunk]              [vir_subtype] [          position               ]

*/


//---------------------------------------------------------------

typedef	struct HASH_struct  HASH;
struct HASH_struct
{
	unsigned int fromPos;
	unsigned int toPos;
};

typedef	struct HEADER_struct  HEADER;
struct HEADER_struct
{
	unsigned int version;
	unsigned int tagcnt;
	unsigned int hashSize;
	unsigned int taglen;
	unsigned long maxUID;
	unsigned long virosaurusDBsize;
	unsigned int maxVirusID;
	unsigned int filterVirusID;
	unsigned int reserved[118];
};

typedef	struct RSLT_struct  RSLT;
struct RSLT_struct
{
	int tag1id;
	int tag2id;
	int tag1pos;
	int tag2pos;
	unsigned int tagpair;
	unsigned int flags;
};

//---------------------------------------------------------------
// vim: tabstop=2 shiftwidth=2

