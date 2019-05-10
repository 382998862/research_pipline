/* 
   Auxiliary software for the project:
   extract-as - discovery of splice variation by comparison of transcripts
   inferred from RNA-seq data 

   Copyright (C) 2011-2013  Liliana Florea (GPL)
*/

#ifndef LIB_H
#define LIB_H

#define min(x,y) (((x)<=(y)) ? (x):(y))
#define max(x,y) (((x)>=(y)) ? (x):(y))

void fatal(const char *msg);
void fatalf(const char *msg, const char *s);
FILE *ckopen(const char *name, const char *mode);
void *ckalloc(size_t amount);
void *ckallocz(size_t amount);
char *strsave(const char *msg);
void *ckrealloc(void * p, size_t size);

#endif
