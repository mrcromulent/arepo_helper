#include "gadgetSnap.h"
#include "rgadget.h"

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

t_gadgetSnap* gc_init( char* filename ) {
  char fn[255];
  int header, footer;
  size_t filesize;
  struct t_gadgetSnapProperty* prop = NULL;
  
  t_gadgetSnap* snap = (t_gadgetSnap*)malloc( sizeof( t_gadgetSnap ) );
  snap->filename = (char*)malloc( strlen( filename ) + 1 );
  strncpy( snap->filename, filename, strlen( filename ) + 1 );
  
  snap->fp = fopen( snap->filename, "r" );
  if (!snap->fp) {
    sprintf( fn, "%s.0", snap->filename );
    snap->fp = fopen( fn, "r" );
    
    if (!snap->fp) {
      printf( "Neither `%s` nor `%s.0` exist.\n", snap->filename, snap->filename );
      return 0;
    }
    
    snap->multipleFiles = 1;
  } else {
    snap->multipleFiles = 0;
  }

  snap->currentFile = 0;
  snap->currentParticle = 0;
  snap->props = 0;
  snap->particlesLoadAtOnce = 0;
  snap->particlesLoaded = 0;
  
  fseek( snap->fp, 0, SEEK_END );
  filesize = ftell( snap->fp );
  fseek( snap->fp, 0, SEEK_SET );
  while ((size_t)ftell( snap->fp ) < filesize) {
    /* add new element at the end of the list */
    if (!snap->props) {
      snap->props = (struct t_gadgetSnapProperty*)malloc( sizeof( struct t_gadgetSnapProperty ) );
      prop = snap->props;
    } else {
      prop->next = (struct t_gadgetSnapProperty*)malloc( sizeof( struct t_gadgetSnapProperty ) );
      prop = prop->next;
    }
    
    prop->next = 0;
    fread( &header, 4, 1, snap->fp );
    fread( prop->name, 4, 1, snap->fp );
    fread( &prop->size, 4, 1, snap->fp );
    fread( &footer, 4, 1, snap->fp );
    
    prop->name[4] = 0;
    prop->offset = ftell( snap->fp ) + 4; /* header */
    prop->size -= 8; /* header + footer */
    prop->dim = 0;
    prop->used = 0;
    prop->datasize = 0;
    prop->data = 0;
    
    fseek( snap->fp, prop->size + 8, SEEK_CUR );
    
    /*printf( "Block %4s of size %d at offset %ld found.\n", prop->name, prop->size, prop->offset );*/
  }
  
  gc_enableprop( snap, "HEAD" );
  gc_readprop( snap, "HEAD", 0, 0 );

  /* reset to read dimensions and calculate some things */
  gc_reset( snap );
  
  return snap;
}

void gc_deinit( t_gadgetSnap* snap ) {
  struct t_gadgetSnapProperty *prop, *next;
  
  if (snap->props) {
    prop = snap->props;
    while (prop->next) {
      if (prop->data)
        free( prop->data );
      
      next = prop->next;
      free( prop );
      
      prop = next;
    }
  }
  
  free( snap );
}

void gc_calcOffsets( t_gadgetSnap* snap ) {
  struct t_gadgetSnapProperty *prop;
  char name[5];
  int size, header, footer;
  size_t filesize;

  fseek( snap->fp, 0, SEEK_END );
  filesize = ftell( snap->fp );
  fseek( snap->fp, 0, SEEK_SET );

  while ((size_t)ftell( snap->fp ) < filesize) {
    fread( &header, 4, 1, snap->fp );
    fread( name, 4, 1, snap->fp );
    fread( &size, 4, 1, snap->fp );
    fread( &footer, 4, 1, snap->fp );
    name[4] = 0;

    prop = snap->props;
    while (prop) {
      if (!strncmp( prop->name, name, 4 )) {
        break;
      }
      prop = prop->next;
    }
  
    if (prop) {
      prop->offset = ftell( snap->fp ) + 4; /* header */
      prop->size = size - 8; /* header + footer */
    }
    
    fseek( snap->fp, size, SEEK_CUR );
  }
}

void gc_reset( t_gadgetSnap* snap ) {
  struct t_gadgetSnapProperty* prop;
  int dims;
  char fn[255];
  
  dims = 0;
  prop = snap->props; 
  while (prop) {
    if (!prop->dim)
      prop->dim = prop->size / ( snap->nparttot * 4 ); /* data is either an array of integer or float */
    dims += prop->dim;
    
    prop = prop->next;
  }
  
  snap->particlesLoadAtOnce = min( snap->npartalltot, MAXDATAMEMORY / ( 4 * dims ) );
  
  if (snap->multipleFiles) {
    if (snap->fp)
      fclose( snap->fp );
    sprintf( fn, "%s.0", snap->filename );
    snap->fp = fopen( fn, "r" );
    
    gc_calcOffsets( snap );
  }
  fseek( snap->fp, 0, SEEK_SET );
  
  snap->currentFile = 0;
  snap->currentParticle = 0;
  snap->currentLocalParticle = 0;
}

void gc_selectFile( t_gadgetSnap* snap, int filenr ) {
  char fn[255];

  if (!snap->multipleFiles) {
    gc_reset( snap );
    return;
  }

  if (filenr >= snap->multipleFiles) {
    printf( "Snapshot does not contain file %d.\n", filenr );
    return;
  }

  if (snap->fp)
    fclose( snap->fp );

  snap->currentFile = filenr;
  sprintf( fn, "%s.%d", snap->filename, snap->currentFile );
  snap->fp = fopen( fn, "r" );
  
  snap->currentLocalParticle = 0;
  
  gc_readprop( snap, "HEAD", 0, 0 );

  gc_calcOffsets( snap );
}

int gc_readnext( t_gadgetSnap* snap ) {
  struct t_gadgetSnapProperty* prop;
  int npart;
  char fn[255];
  
  if (snap->currentParticle == snap->npartalltot)
    return 0;
  
  npart = min( snap->particlesLoadAtOnce, snap->npartalltot - snap->currentParticle );
  
  if (snap->multipleFiles && snap->nparttot > snap->currentLocalParticle + npart) {
    if (snap->currentLocalParticle == snap->nparttot) {
      if (snap->fp)
        fclose( snap->fp );
      snap->currentFile++;
      sprintf( fn, "%s.%d", snap->filename, snap->currentFile );
      snap->fp = fopen( fn, "r" );
      
      snap->currentLocalParticle = 0;
      
      gc_readprop( snap, "HEAD", 0, 0 );
      
      gc_calcOffsets( snap );
      
      npart = min( npart, snap->nparttot );
    } else {
      npart = snap->nparttot - snap->currentLocalParticle;
    }
  }
  
  prop = snap->props;
  while (prop) {
    if (prop->used)
      gc_readprop( snap, prop->name, snap->currentLocalParticle, npart );
    prop = prop->next;
  }
  
  snap->currentLocalParticle += npart;
  snap->currentParticle += npart;
  snap->particlesLoaded = npart;
  
  return npart;
}

int gc_readprop( t_gadgetSnap* snap, char *propname, int startparticle, int nparticles ) {
  struct t_gadgetSnapProperty* prop;
  int i;
  
  prop = snap->props;
  while (prop) {
    if (!strcmp( prop->name, propname )) {
      break;
    }
    prop = prop->next;
  }
  
  if (!prop) return 0; /* property does not exist */
  
  if (!strcmp( propname, "HEAD" )) {
    fseek( snap->fp, prop->offset, SEEK_SET );
    fread( snap->npart, 4, 6, snap->fp );
    fread( snap->masses, 8, 6, snap->fp );
    fread( &snap->time, 8, 1, snap->fp );
    fread( &snap->redshift, 8, 1, snap->fp );
    fread( &snap->flag_sfr, 4, 1, snap->fp );
    fread( &snap->flag_feedback, 4, 1, snap->fp );
    fread( snap->npartall, 4, 6, snap->fp );
    fread( &snap->flag_cooling, 4, 1, snap->fp );
    fread( &snap->multipleFiles, 4, 1, snap->fp );
    fseek( snap->fp, 128, SEEK_CUR ); /* skip 128 (empty) bytes */
    
    snap->nparttot = 0;
    snap->npartalltot = 0;
    for (i=0; i<6; i++) {
      snap->nparttot += snap->npart[i];
      snap->npartalltot += snap->npartall[i];
    }
    
    if (!snap->npartalltot)
      snap->npartalltot = snap->nparttot; /* npartalltot sometimes not set for single file output */

    if (snap->multipleFiles == 1)
                  snap->multipleFiles = 0;                /* only > 0 if there is more than one file */
  } else {
    if (!prop->dim)
      prop->dim = prop->size / ( snap->nparttot * 4 ); /* data is either an array of integer or float */
    
    if (prop->datasize != nparticles) {
      if (prop->data)
        free( prop->data );
      prop->datasize = nparticles;
      prop->data = (char*)malloc( prop->datasize * prop->dim * 4 );

      if (!prop->data)
              printf( "Could not allocate %d bytes for %s\n", prop->datasize * prop->dim * 4, prop->name );
    }
    
    fseek( snap->fp, prop->offset + startparticle * prop->dim * 4, SEEK_SET );
    fread( prop->data, 4, nparticles * prop->dim, snap->fp );
  }
  
  return 1;
}

void gc_freeprop( t_gadgetSnap* snap, char *propname ) {
  struct t_gadgetSnapProperty* prop;
  
  prop = snap->props;
  while (prop) {
    if (!strcmp( prop->name, propname )) {
      break;
    }
    prop = prop->next;
  }
  
  if (!prop) return; /* property does not exist */

  if (prop->data) {
    free( prop->data );
    prop->data = 0;
    prop->datasize = 0;
  }
}

char* gc_getdata( t_gadgetSnap* snap, char *propname ) {
  struct t_gadgetSnapProperty* prop;
  prop = snap->props;
  
  while (prop) {
    if (!strcmp( prop->name, propname )) {
      return prop->data;
    }
    prop = prop->next;
  }
  
  return 0;
}

int gc_getdim( t_gadgetSnap* snap, char *propname ) {
  struct t_gadgetSnapProperty* prop;
  prop = snap->props;
  
  while (prop) {
    if (!strcmp( prop->name, propname )) {
      return prop->dim;
    }
    prop = prop->next;
  }
  
  return 0;
}

int gc_getpropcount( t_gadgetSnap* snap ) {
  int count;
  struct t_gadgetSnapProperty* prop;
  prop = snap->props;
  
  count = 0;
  while (prop) {
    count++;
    prop = prop->next;
  }
  
  return count;
}

int gc_getpropname( t_gadgetSnap* snap, int propid, char* propname ) {
  int count;
  struct t_gadgetSnapProperty* prop;
  prop = snap->props;
  
  count = 0;
  while (prop) {
    if (propid == count) {
      memcpy( propname, prop->name, 5 );
      return 0;
    }
    count++;
    prop = prop->next;
  }
  
  return -1;
}

int gc_hasprop( t_gadgetSnap* snap, char *propname ) {
  struct t_gadgetSnapProperty* prop;
  prop = snap->props;
  
  while (prop) {
    if (!strcmp( prop->name, propname )) {
      return 1;
    }
    prop = prop->next;
  }
  
  return 0;
}

void gc_disableprop( t_gadgetSnap* snap, char *propname ) {
  struct t_gadgetSnapProperty* prop;
  
  if (!strcmp( "HEAD", propname)) {
    printf( "`HEAD` must always be enabled.\n" );
    return;
  }
  
  prop = snap->props;
  while (prop) {
    if (!strcmp( prop->name, propname )) {
      prop->used = 0;
      return;
    }
    prop = prop->next;
  }
}

void gc_enableprop( t_gadgetSnap* snap, char *propname ) {
  struct t_gadgetSnapProperty* prop;

  prop = snap->props; 
  while (prop) {
    if (!strcmp( prop->name, propname )) {
      prop->used = 1;
      return;
    }
    prop = prop->next;
  }
}
