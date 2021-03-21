#ifndef GADGETSNAP_H
#define GADGETSNAP_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct t_gadgetSnapProperty {
	char name[5];
	size_t offset;
	int size;
	int dim;
	int used;
	int datasize;
	char *data;
	struct t_gadgetSnapProperty* next;
};

typedef struct {
	char *filename;
	FILE* fp;
	int multipleFiles;
	int currentFile;
	int currentParticle;
	int currentLocalParticle;
	int particlesLoadAtOnce;
	int particlesLoaded;
	struct t_gadgetSnapProperty *props;
	
	double time, redshift;
	int npart[6], nparttot;		/* particles in current file */
	double masses[6];
	int flag_sfr, flag_feedback, flag_cooling;
	int npartall[6], npartalltot;	/* particles in all files of this snapshot */
} t_gadgetSnap;

t_gadgetSnap* gc_init( char* filename );
void gc_deinit( t_gadgetSnap* snap );
void gc_calcOffsets( t_gadgetSnap* snap );
void gc_reset( t_gadgetSnap* snap );
void gc_selectFile( t_gadgetSnap* snap, int filenr );
int gc_readnext( t_gadgetSnap* snap );
int gc_readprop( t_gadgetSnap* snap, char *propname, int startparticle, int nparticles );
void gc_freeprop( t_gadgetSnap* snap, char *propname );
char* gc_getdata( t_gadgetSnap* snap, char *propname );
int gc_getdim( t_gadgetSnap* snap, char *propname );
int gc_getpropcount( t_gadgetSnap* snap );
int gc_getpropname( t_gadgetSnap* snap, int propid, char* propname );
int gc_hasprop( t_gadgetSnap* snap, char *propname );
void gc_disableprop( t_gadgetSnap* snap, char *propname );
void gc_enableprop( t_gadgetSnap* snap, char *propname );

#endif
