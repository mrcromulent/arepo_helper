#ifndef RGADGET_H
#define RGADGET_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXDATAMEMORY  1024000000	/* 1 GB */
#define MAXPARTPERICFILE 20000000	/* 2e7  */

void gadget_writeHeader( FILE* fd, int npartall, int npart, int num_files );
void gadget_writeBlockFloat( FILE* fd, char *name, int npart, int dim, float *data );
void gadget_writeBlockInt( FILE* fd, char *name, int npart, int dim, int *data );


#endif /* RGADGET_H */
