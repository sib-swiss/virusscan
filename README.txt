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
