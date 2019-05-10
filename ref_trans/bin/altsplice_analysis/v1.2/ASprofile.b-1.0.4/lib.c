/* 
   Auxiliary software for the project:
   extract-as - discovery of splice variation by comparison of transcripts
   inferred from RNA-seq data 

   Copyright (C) 2011-2013  Liliana Florea (GPL)
*/

#include "libc.h"
#include "lib.h"

/* Auxiliary functions */

/* fatal ---------------------------------------------- print message and die */
void fatal(const char *msg)
{
     fprintf(stderr, "%s\n", msg);
     exit(1);
}


/* fatalf ---------------------------------------------- print message and die */
void fatalf(const char *msg, const char *s)
{
     char *buf = (char *)ckalloc(strlen(msg)+strlen(s)+1);
     sprintf(buf, msg, s);
     fatal(buf);
}

/* ckopen -------------------------------------- open file; check for success */
FILE *ckopen(const char *name, const char *mode)
{
        FILE *fp;

        if ((fp = fopen(name, mode)) == NULL)
                fatalf("Cannot open %s.", name);
        return fp;
}

/* ckalloc -------------------------------- allocate space; check for success */
void *ckalloc(size_t amount)
{
        void *p;
        char buf[100];

        if ((long)amount<0) 
            printf("Amount %ld %ld\n", amount, (long)amount);
        assert((long)amount >= 0);
        if (amount == 0)
                amount = 1; /* ANSI portability hack */
        if ((p = malloc(amount)) == NULL) {
             sprintf(buf, "Ran out of memory trying to allocate %lu.",
                        (unsigned long)amount);
             fatal(buf);
        }
#if 0
        memset(p, 0, amount); /* XXX */
#endif 
        return p;
}

/* ckallocz -------------------- allocate space; zero fill; check for success */
void *ckallocz(size_t amount)
{
        void *p = ckalloc(amount);
        memset(p, 0, amount);
        return p;
}

/* strsave ----------------------- allocate space for msg and store msg there */ 
char *strsave(const char *msg)
{
   int amount = strlen(msg);
   char *p = (char *)ckalloc(amount+1); 
   strcpy(p,msg); p[amount] = '\0';

   return p;
}

/* ckrealloc ---------------------- check and reallocate space as needed */
void *ckrealloc(void * p, size_t size)
{
        p = p ? realloc(p, size) : malloc(size);
        if (!p)
                fatal("ckrealloc failed");
        return p;
}
