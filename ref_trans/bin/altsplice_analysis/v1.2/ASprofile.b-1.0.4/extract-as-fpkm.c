#include "libc.h"
#include "lib.h"

#define DEFAULT_ALLOC        10000
#define DEFAULT_ALLOC_SMALL    100

#define DEFAULT_W                0

/* Compiling: gcc extract_as_fpkm.c lib.c -o extract_as_fpkm -lm */

#define V      100

static char Usage[] = "%s txptgtf genome_hdrs eventsnr [-W wiggle]\n";

static char **Chroms;
static int    numChroms = 0;
static int    maxChroms = 0;

typedef struct {
   char  transcript_id[100];
   char  gene_id[100];
   char  gene_name[100];
   int   gaid;
   char  ori;
   int   num_exons;
   int   max_exons;
   int  *from, *to;
   char *text;
   double fpkm;
   int   printed;
   int   len;
} transcript_t;


#define TSS_TYPE            11
#define TTS_TYPE            12
#define SKIP_TYPE           20
#define SKIP_ON_TYPE        21
#define SKIP_OFF_TYPE       22
#define MSKIP_ON_TYPE       23
#define MSKIP_OFF_TYPE      24
#define XSKIP_ON_TYPE       25
#define XSKIP_OFF_TYPE      26
#define XMSKIP_ON_TYPE      27
#define XMSKIP_OFF_TYPE     28
#define IR_TYPE             30
#define IR_ON_TYPE          31
#define IR_OFF_TYPE         32
#define MIR_ON_TYPE         33
#define MIR_OFF_TYPE        34
#define XIR_ON_TYPE         35
#define XIR_OFF_TYPE        36
#define XMIR_ON_TYPE        37
#define XMIR_OFF_TYPE       38
#define AE_TYPE             40
#define XAE_TYPE            41

static char *Codes[50] = { "", "", "", "", "", "", "", "", "", "",
                           "", "TSS", "TTS", "", "", "", "", "", "", "",
                           "SKIP", "SKIP_ON", "SKIP_OFF", "MSKIP_ON", "MSKIP_OFF",
                           "XSKIP_ON", "XSKIP_OFF", "XMSKIP_ON", "XMSKIP_OFF", "",
                           "IR", "IR_ON", "IR_OFF", "MIR_ON", "MIR_OFF",
                           "XIR_ON", "XIR_OFF", "XMIR_ON", "XMIR_OFF", "",
                           "AE", "XAE", "", "", "", "", "", "", "", ""};

static int numCodes = 50;

typedef struct event {
  int  event_id;
  int  type;
  int  gaid, from, to;
  char pattern[10000];
  char ori;
  char gene_id[100];
  double fpkm;
  char other[100];
} event_t;

typedef struct junc {
  int  from, to, gaid, num;
  char ori;
} junc_t;

static int  lookup_chrom(char *,char *chroms[],int);
static int  lookup_event_type(char *);
static int  junc_lookup2(int,junc_t *,int,int);
static void load_chroms(char *gafile, char ***chroms, int *numchroms, int *maxchroms);
#if 0
static void loadTranscripts(const char *,transcript_t **,int *,int *);
#endif
static void loadGTF(const char *,transcript_t **,int *,int *);
static void loadASVs(const char *,transcript_t **,int *, int *);
static void loadJuncs(const char *,junc_t **,int *, int*);
static int  transcript_cmp(const void *,const void *);
static int  junc_cmp(const void *,const void *);
static int  overlap(int,int,int,int);
static int  readEventLine(char *,event_t *);
static int  searchIntronApprox(int,int,transcript_t *,int,int *);
static int  searchExonApprox(int,int,transcript_t *,int,int *);
static int  searchJuncApprox(int,int,int,junc_t *,int);


static void transcript2text(transcript_t *tptr, char text[])
{
   char  *s = text;
   int    st, j;

   st = sprintf(s, "%d-%d", tptr->from[0], tptr->to[0]);
   s += st;

   for (j=1; j<tptr->num_exons; j++) {
     st = sprintf(s, ",%d-%d", tptr->from[j], tptr->to[j]);
     s += st;
   }
   *s = '\0';
}

void main(int argc, char *argv[])
{
    transcript_t **Transcripts, *tptr;
    int           *transcriptCounts;
    int            numTranscripts, maxTranscripts;
    junc_t       **Juncs;
    int           *juncCounts = NULL;
    int            numJuncs, maxJuncs;
    int            i, t, W, arg, MaxTranscriptSpan, isASV;
    event_t        ev;
    FILE          *fp;
    char           buf[100000], *tok, juncFile[100];
  
    if (argc<4) fatalf(Usage, argv[0]);

    W = DEFAULT_W;

    arg = 4;
    isASV = 0;
    *juncFile = '\0';
    while (arg < argc) {
      if (!strncmp(argv[arg], "-asv", 4))
        isASV = 1; 
      else if (!strncmp(argv[arg], "-W", 2))
        W = atoi(argv[++arg]);
      else if (!strncmp(argv[arg], "-j",2))
        strcpy(juncFile, argv[++arg]);
      else 
        fatalf(Usage, argv[0]);
      arg++;
    }

    load_chroms(argv[2], &Chroms, &numChroms, &maxChroms);

    numTranscripts = maxTranscripts = 0;

    Transcripts = (transcript_t **)malloc(numChroms*sizeof(transcript_t));
    memset(Transcripts, 0, numChroms*sizeof(transcript_t *));

    if (!isASV) 
      loadGTF(argv[1], &(Transcripts[0]), &numTranscripts, &maxTranscripts);
    else 
      loadASVs(argv[1], &(Transcripts[0]), &numTranscripts, &maxTranscripts);

    qsort(Transcripts[0], numTranscripts, sizeof(transcript_t), transcript_cmp);
    transcriptCounts = (int *)malloc(numChroms*sizeof(int));
    memset(transcriptCounts, 0, numChroms*sizeof(int));

    for (i=0; i<numTranscripts; i++) {
      int gaid = Transcripts[0][i].gaid;
      if (Transcripts[gaid]==NULL) Transcripts[gaid] = Transcripts[0]+i;
      transcriptCounts[gaid]++;
    }

    for (i=0, tptr=Transcripts[0], MaxTranscriptSpan=0; i<numTranscripts; i++)
      if ((t=tptr->to[tptr->num_exons-1]-tptr->from[0]+1) > MaxTranscriptSpan)
        MaxTranscriptSpan = t; 


    Juncs = NULL;
    numJuncs = maxJuncs = 0;

    if (*juncFile) {
      Juncs = (junc_t **)malloc(numChroms*sizeof(junc_t));
      memset(Juncs, 0, numChroms*sizeof(junc_t *));

      loadJuncs(juncFile, &(Juncs[0]), &numJuncs, &maxJuncs);

      qsort(Juncs[0], numJuncs, sizeof(junc_t), junc_cmp);
      juncCounts = (int *)malloc(numChroms*sizeof(int));
      memset(juncCounts, 0, numChroms*sizeof(int));

      for (i=0; i<numJuncs; i++) {
        int gaid = Juncs[0][i].gaid;
        if (Juncs[gaid]==NULL) Juncs[gaid] = Juncs[0]+i;
        juncCounts[gaid]++;
      }
    }

    /* print header */
    printf("event_id\tevent_type\tgene_id\tchrom\tevent_start\tevent_end\tevent_pattern\tstrand\tfpkm\n");

    fp = ckopen(argv[3], "r");
    while (fgets(buf, 100000, fp)!=NULL) {

      if (!strncmp(buf,"event_id",8)) continue;

      int st = readEventLine(buf, &ev);
      if (st) {
        fprintf(stderr, "Warning: Something wrong with event id %s. Skipping.\n", buf);
        continue; 
      }

      if (Juncs!=NULL) {
        int lf, tr, lr, l, f, t, r;
        switch (ev.type) {
          case SKIP_ON_TYPE:
               if (sscanf(ev.pattern, "%d,%d-%d,%d", &l, &f, &t, &r)!=4)
                 fatalf("Error: Something wrong with SKIP_ON event pattern. %s\n", ev.pattern);
               lf = searchJuncApprox(l,f,W,Juncs[ev.gaid],juncCounts[ev.gaid]);
               tr = searchJuncApprox(t,r,W,Juncs[ev.gaid],juncCounts[ev.gaid]);
               lr = searchJuncApprox(l,r,W,Juncs[ev.gaid],juncCounts[ev.gaid]);
               sprintf(ev.other, "%d,%d,%d", lf, tr, lr);

               break;

          case SKIP_OFF_TYPE:
               if (sscanf(ev.pattern, "%d,%d", &l, &r)!=2)
                 fatalf("Error: Something wrong with SKIP_OFF event pattern. %s\n", ev.pattern);
               lr = searchJuncApprox(l,r,W,Juncs[ev.gaid],juncCounts[ev.gaid]);
               sprintf(ev.other, "%d", lr);

               break;

          default:
               break;
        }
      }

      for (i=0, tptr=Transcripts[ev.gaid]; i<transcriptCounts[ev.gaid]; i++, tptr++) {

        if (tptr->to[tptr->num_exons-1]<ev.from) continue;
        if (ev.to<tptr->from[0]) break;

        if (overlap(ev.from,ev.to,tptr->from[0],tptr->to[tptr->num_exons-1])) {
          switch (ev.type) {
            case TSS_TYPE: if (tptr->num_exons>1) {
                             if ((tptr->ori=='+' && abs(tptr->to[0]-atoi(ev.pattern))<=W) ||
                                 (tptr->ori=='-' && abs(tptr->from[tptr->num_exons-1]-atoi(ev.pattern))<=W))
                               ev.fpkm += tptr->fpkm;
                           }
                           break;
                            
            case TTS_TYPE: if (tptr->num_exons>1) {
                             if ((tptr->ori=='-' && abs(tptr->to[0]-atoi(ev.pattern))<=W) ||
                                 (tptr->ori=='+' && abs(tptr->from[tptr->num_exons-1]-atoi(ev.pattern))<=W))
                               ev.fpkm += tptr->fpkm;
                           }
                           break;

            case SKIP_ON_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int l, f, t, r, xi, st = 0;
                             if (sscanf(ev.pattern, "%d,%d-%d,%d", &l, &f, &t, &r)!=4)
                               fatalf("Error: Something wrong with SKIP_ON event pattern. %s\n", ev.pattern);
                             st = searchExonApprox(f,t,tptr,W,&xi);
                             if (!st && 
                                 (xi>=1) && abs(tptr->to[xi-1]-l)<=W &&
                                 (xi<tptr->num_exons-1) && abs(tptr->from[xi+1]-r)<=W)
                               ev.fpkm += tptr->fpkm;
                           }

                           break;

            case SKIP_OFF_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int f, t, yi, st = 0;
                             if (sscanf(ev.pattern, "%d,%d", &f, &t)!=2)
                               fatalf("Error: Something wrong with SKIP_OFF event pattern. %s\n", ev.pattern);
                             st = searchIntronApprox(f,t,tptr,W,&yi);
                             if (!st) ev.fpkm += tptr->fpkm;
                           }
                           break;

            case MSKIP_ON_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int l, f, t, r, yi, st = 0;
                             char *s = ev.pattern;
                             if (sscanf(ev.pattern, "%d,%d-%d,%d", &l, &f, &t, &r)!=4)
                               fatalf("Error: Something wrong with MSKIP_ON event pattern (1). %s\n", ev.pattern);
                             st = searchIntronApprox(l,f,tptr,W,&yi);
                             while (*s && isdigit(*s)) s++; s++;
                             while (*s && isdigit(*s)) s++;

                             while (!st && *s && abs(tptr->to[yi]-t)<=W && 
                                    (yi+1<tptr->num_exons) && abs(tptr->from[yi+1]-r<=W)) {
                                s++; while (*s && isdigit(*s)) s++;
                                s++; while (*s && isdigit(*s)) s++;
                               
                                if (*s=='-') {
                                  if (sscanf(s, "-%d,%d", &t, &r)!=2)
                                    fatalf("Error: Something wrong with MSKIP_ON event pattern (2). %s\n", s);
                                  yi++; 
                                } else
                                  ev.fpkm += tptr->fpkm;
                              }
                           }
                           break;

            case MSKIP_OFF_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int f, t, yi, st = 0;
                             if (sscanf(ev.pattern, "%d,%d", &f, &t)!=2)
                               fatalf("Error: Something wrong with MSKIP_OFF event pattern. %s\n", ev.pattern);
                             st = searchIntronApprox(f,t,tptr,W,&yi);
                             if (!st) ev.fpkm += tptr->fpkm;
                           }
                           break;

            case IR_ON_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int f, t, xi, st = 0;
                             if (sscanf(ev.pattern, "%d-%d", &f, &t)!=2)
                               fatalf("Error: Something wrong with IR_ON event pattern. %s\n", ev.pattern);
                             st = searchExonApprox(f,t,tptr,W,&xi);
                             if (!st && xi>0 && xi<tptr->num_exons-1) ev.fpkm += tptr->fpkm;
                           }
                           break;

            case IR_OFF_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int lb, le, rb, re, xi, st = 0;
                             if (sscanf(ev.pattern, "%d-%d,%d-%d", &lb, &le, &rb, &re)!=4)
                               fatalf("Error: Something wrong with IR_OFF event pattern. %s\n", ev.pattern);
                             st = searchExonApprox(lb,le,tptr,W,&xi);
                             if (!st && (xi<tptr->num_exons-1) &&
                                 abs(tptr->from[xi+1]-rb)<=W && abs(tptr->to[xi+1]-re)<=W)
                               ev.fpkm += tptr->fpkm;
                           }
                           break;

            case MIR_ON_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int f, t, xi, st = 0;
                             if (sscanf(ev.pattern, "%d-%d", &f, &t)!=2)
                               fatalf("Error: Something wrong with MIR_ON event pattern. %s\n", ev.pattern);
                             st = searchExonApprox(f,t,tptr,W,&xi);
                             if (!st && xi>0 && xi<tptr->num_exons-1) ev.fpkm += tptr->fpkm;
                           }
                           break;

            case MIR_OFF_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int lb, le, rb, re, xi, st = 0;
                             char *s = ev.pattern;
                             if (sscanf(ev.pattern, "%d-%d,%d-%d", &lb, &le, &rb, &re)!=4)
                               fatalf("Error: Something wrong with MIR_OFF event pattern (1). %s\n", ev.pattern);
                             st = searchExonApprox(lb,le,tptr,W,&xi);
                             while (*s && isdigit(*s)) s++; s++;
                             while (*s && isdigit(*s)) s++;

                             while (!st && *s && (xi<tptr->num_exons-1) &&
                                    abs(tptr->from[xi+1]-rb)<=W && abs(tptr->to[xi+1]-re)<=W) {
                               s++; while (*s && isdigit(*s)) s++;
                               s++; while (*s && isdigit(*s)) s++;

                               if (*s==',') {
                                 if (sscanf(s, ",%d-%d", &rb, &re)!=2)
                                   fatalf("Error: Something wrong with MIR_OFF event pattern (2). %s\n", s);
                                 xi++;
                               } else
                                 ev.fpkm += tptr->fpkm;
                             }
                           }
                           break;

            case AE_TYPE:
                           tok = strstr(tptr->text, ev.pattern);
                           if (!W && tok) ev.fpkm += tptr->fpkm;
                           else {
                             int f, t, xi, st = 0;
                             if (sscanf(ev.pattern, "%d-%d", &f, &t)!=2)
                               fatalf("Error: Something wrong with AE event pattern. %s\n", ev.pattern);
                             st = searchExonApprox(f,t,tptr,W,&xi);
                             if (!st && xi>0 && xi<tptr->num_exons-1) ev.fpkm += tptr->fpkm;
                           }
                           break;

            default:       tok = strstr(tptr->text, ev.pattern);
                           if (tok) ev.fpkm += tptr->fpkm; 
                           break;
          } 
        }
      }
      printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%1.10f",
             ev.event_id, Codes[ev.type], ev.gene_id, Chroms[ev.gaid],
             ev.from, ev.to, ev.pattern, ev.ori, ev.fpkm);

      if (*juncFile) printf("\t%s\n", ev.other);
      else printf("\n"); 
    }
    fclose(fp);

    for (i=0; i<numTranscripts; i++) {
      free(Transcripts[0][i].text);
      free(Transcripts[0][i].from);
      free(Transcripts[0][i].to);
    }
    free(Transcripts[0]);
    free(Transcripts);

    if (Juncs) {
      free(Juncs[0]); free (Juncs);
    }

    free(Chroms);
}


static void load_chroms(char *gafile, char ***chroms, int *numchroms, int *maxchroms)
{
  char  buf[10000], *tok;
  FILE *fp = ckopen(gafile, "r");

  *maxchroms = *numchroms = 0;
  *chroms = NULL;

  while (fgets(buf,10000,fp)!=NULL) {
    if (*buf!='>') continue;

    if (*numchroms>=*maxchroms) {
      *maxchroms = max(2*(*maxchroms),DEFAULT_ALLOC_SMALL);
      *chroms = realloc(*chroms,(*maxchroms)*sizeof(char *));
    }

    tok = strtok(buf, " \t\n");
    (*chroms)[*numchroms] = strsave(tok+1);

    (*numchroms)++;
  }
  fclose(fp);
}

static void loadASVs(const char *filename, transcript_t **transcripts, int *numtranscripts, int *maxtranscripts)
{
   FILE         *fp = ckopen(filename, "r");
   char          buf[500000];
   int           from, to, idummy;
   char         *tok;
   transcript_t *txpt;

   *numtranscripts = *maxtranscripts = 0;

   next_asv:
   while (fgets(buf, 100000, fp)!=NULL) {
     /* ASVid gaid geneid strand ori exon_list */

     if (*numtranscripts>=*maxtranscripts) {
       *maxtranscripts = max(DEFAULT_ALLOC, 2*(*maxtranscripts));
       *transcripts = realloc(*transcripts, (*maxtranscripts)*sizeof(transcript_t));
     }

     txpt = (*transcripts)+(*numtranscripts);

     tok = strtok(buf, " ");
     strcpy(txpt->transcript_id, tok);
     tok = strtok(NULL, " ");
     txpt->gaid = atoi(tok); 
     tok = strtok(NULL, " ");
     strcpy(txpt->gene_id, tok);
     strcpy(txpt->gene_name, tok);
     tok = strtok(NULL, " ");
     if (!strcmp(tok, "-1"))
       txpt->ori = '-';
     else if (!strcmp(tok, "1"))
       txpt->ori = '+';
     else
       txpt->ori = '.';
     tok = strtok(NULL, " ");  /* skip */

     txpt->fpkm = 1;

     /* now process the exons */
     txpt->num_exons = txpt->max_exons = txpt->len = 0;
     txpt->from = txpt->to = NULL;
     txpt->printed = 0;

     tok = strtok(NULL, " \n");
     while (tok) {
       if (txpt->num_exons>=txpt->max_exons) {
         txpt->max_exons = max(DEFAULT_ALLOC_SMALL, 2*txpt->max_exons);
         txpt->from = realloc(txpt->from, txpt->max_exons*sizeof(int));
         txpt->to = realloc(txpt->to, txpt->max_exons*sizeof(int));
       }

       if (sscanf(tok, "ND_%d(%d-%d)", &idummy, &from, &to)!=3) {
         fprintf(stderr, "Warning: Unrecognized ASV line. Skipping ASV. %s\n", txpt->transcript_id);
         goto next_asv; 
       }

       txpt->from[txpt->num_exons] = from;
       txpt->to[txpt->num_exons] = to;

       txpt->num_exons++;

       tok = strtok(NULL, " \n");
     }
    
     assert(txpt->num_exons!=0);
     txpt->text = ckalloc(txpt->num_exons*(9*2+2)*sizeof(char));
     transcript2text(txpt, txpt->text);

     (*numtranscripts)++;
   }
   fclose(fp);
}

static void loadGTF(const char *filename, transcript_t **goldrefs, int *numgoldrefs, int *maxgoldrefs)
{
   FILE         *fp = ckopen(filename, "r");
   char          buf[100000], chrom[20], ori;
   int           i, from, to, idxc;
   char         *tok, sdummy1[100], sdummy2[100], sdummy3[100], prev_transcript_id[100], info[100000];
   double        fpkm = -1.0;
   transcript_t *txpt;

   *numgoldrefs = *maxgoldrefs = 0;

   while (fgets(buf, 100000, fp)!=NULL) {

      if (buf[0] == '#') continue;

      /* chr3  Cufflinks       transcript      71002043        71114061        1000    -       .       gene_id "ENSG00000114861"; transcript_id "adipose.289427.1"; FPKM "19.8831"; frac "1.000000"; conf_lo "7.135536"; conf_hi "22.544681"; cov "81.315579";
      */

      /* NT_166433   protein_coding  exon    11955   12166   .    +    .    gene_id "ENSMUSG00000000702"; transcript_id "ENSMUST00000105216"; exon_number "1"; gene_name "AC007307.1"; transcript_name "AC007307.1-201";
      */

      if (sscanf(buf, "%s\t%s\t%s\t%d\t%d\t%s\t%c\t.\t",
                       chrom, sdummy1, sdummy2, &from, &to, sdummy3, &ori) != 7) {
        fprintf(stderr, "Warning: Irregular line in reference file. Skipping line. %s\n", buf);
        continue;
      }
      if (strcmp(sdummy2,"exon"))
        continue;

      if ((idxc=lookup_chrom(chrom, Chroms, numChroms))<0) {
        fprintf(stderr, "Warning: cannot find chromosome. Skipping line. %s\n", buf);
        continue;
      }

      strcpy(info, strstr(buf, "transcript_id")+strlen("transcript_id \""));
      tok = strtok(info, " \n");
      tok[strlen(tok)-2] = '\0';
      strcpy(sdummy2, tok);

      if (strcmp(prev_transcript_id, sdummy2)) {
        /* initialize a new goldref transcript */

        (*numgoldrefs)++;

        if (*numgoldrefs>=*maxgoldrefs) {
          *maxgoldrefs = max(DEFAULT_ALLOC, 2*(*maxgoldrefs));
          *goldrefs = realloc(*goldrefs, (*maxgoldrefs)*sizeof(transcript_t));
        }
        txpt = (*goldrefs)+(*numgoldrefs-1);

        strcpy(info, strstr(buf, "gene_id"));

        *sdummy1 = *sdummy3 = '\0';
        tok = strtok(info, " \n");
        while (tok) {
          if (!strcmp(tok, "gene_id")) {
            tok = strtok(NULL, " \n");
            strncpy(sdummy1, tok+1, strlen(tok)-1);
            sdummy1[strlen(tok)-3] = '\0';
          } else if (!strcmp(tok, "gene_name")) {
            tok = strtok(NULL, " \n");
            strncpy(sdummy3, tok+1, strlen(tok)-1);
            sdummy3[strlen(tok)-3] = '\0';
          } else if (!strcmp(tok, "FPKM") || !strcmp(tok, "fpkm")) {
            tok = strtok(NULL, " \n");
            sscanf(tok+1, "%lf\"", &fpkm);
          }
          tok = strtok(NULL, " \n");
        }

        txpt->gaid = idxc;
        txpt->ori = ori;
        txpt->fpkm = fpkm;
        strcpy(txpt->gene_id, sdummy1);
        strcpy(txpt->transcript_id, sdummy2);
        strcpy(txpt->gene_name, sdummy3);

        txpt->num_exons = txpt->max_exons = 0;
        txpt->from = txpt->to = NULL;
        txpt->len = 0;
        txpt->printed = 0;

        strcpy(prev_transcript_id, sdummy2);

      } else
        txpt = (*goldrefs)+(*numgoldrefs-1);

      /* now fill in the exon */
      if (txpt->num_exons>=txpt->max_exons) {
        txpt->max_exons = max(DEFAULT_ALLOC_SMALL, 2*txpt->max_exons);
        txpt->from = realloc(txpt->from, txpt->max_exons*sizeof(int));
        txpt->to = realloc(txpt->to, txpt->max_exons*sizeof(int));
      }
      txpt->from[txpt->num_exons] = from;
      txpt->to[txpt->num_exons] = to;
      txpt->len += to-from+1;

      txpt->num_exons++;
   }

   fclose(fp);

   /* reorder exons, if necessary (e.g., on the reverse strand) */

   for (i=0; i<*numgoldrefs; i++) {
     txpt = (*goldrefs)+i;

     if ((txpt->num_exons>1) && (txpt->ori=='-') && (txpt->from[0]>txpt->from[2])) {
       int u, l, aux;

       l = 0; u = txpt->num_exons-1;
       while (l<u) {
         aux = txpt->from[l]; txpt->from[l] = txpt->from[u]; txpt->from[u] = aux;
         aux = txpt->to[l]; txpt->to[l] = txpt->to[u]; txpt->to[u] = aux;

         l++; u--;
       }
     }
   }

   /* fill in the text field for all transcripts */
   for (i=0; i<*numgoldrefs; i++) {
     txpt = (*goldrefs)+i;
     txpt->text = ckalloc(txpt->num_exons*(9*2+2)*sizeof(char));
     transcript2text(txpt, txpt->text);
  }

}

#if 0
static void loadTranscripts(const char *filename, transcript_t **transcripts, int *numtranscripts, int *maxtranscripts)
{
   FILE         *fp = ckopen(filename, "r");
   char          buf[100000], chrom[20], ori;
   int           from, to, score, idxc;
   char         *tok, sdummy1[100], sdummy2[100], sdummy3[100], info[100000];
   transcript_t *txpt;
   double        fpkm;

   *numtranscripts = *maxtranscripts = 0;

   while (fgets(buf, 100000, fp)!=NULL) {

      next_transcript:

      if (buf[0] == '#') continue;

      /* chr1    Cufflinks       transcript      4797676 4836758 1000    +       .       gene_id ""; transcript_id "aB.1353.1"; FPKM "7.9807614546"; frac "1.000000"; conf_lo "2.330713"; conf_hi "13.630810"; cov "89.076212";
      */

      if (sscanf(buf, "%s\tCufflinks\ttranscript\t%d\t%d\t%d\t%c\t.\t",
                       chrom, &from, &to, &score, &ori) != 5) {
        fprintf(stderr, "Warning: Unrecognized transcript line. Skipping transcript. %s\n", buf);
        while (fgets(buf, 100000, fp)!=NULL) {
          if (strstr(buf, "\ttranscript\t"))
            goto next_transcript;
        }
        return;  /* reached end of file */
      }
      tok = buf + strlen(chrom) + 22;
      while (isdigit(*tok)) tok++; tok++;
      while (isdigit(*tok)) tok++; tok++;
      while (isdigit(*tok)) tok++; tok+=5;
      strcpy(info, tok);

      /* read in the info part of the line */
      *sdummy1 = *sdummy2 = *sdummy3 = '\0';
      tok = strtok(info, " \n");
      while (tok) {
        if (!strcmp(tok, "gene_id")) {
          tok = strtok(NULL, " \n");
          strncpy(sdummy1, tok+1, strlen(tok)-1);
          sdummy1[strlen(tok)-3] = '\0';
        } else if (!strcmp(tok, "transcript_id")) {
          tok = strtok(NULL, " \n");
          strncpy(sdummy2, tok+1, strlen(tok)-1);
          sdummy2[strlen(tok)-3] = '\0';
        } else if (!strcmp(tok, "gene_name")) {
          tok = strtok(NULL, " \n");
          strncpy(sdummy1, tok+1, strlen(tok)-1);
          sdummy3[strlen(tok)-3] = '\0';
        } else if (!strcmp(tok, "FPKM")) {
          tok = strtok(NULL, " \n");
          sscanf(tok+1, "%lf\";", &fpkm);
        }
        tok = strtok(NULL, " \n");
      }

      if (!*sdummy1 || !*sdummy2)  {
        fprintf(stderr, "Warning: Unrecognized transcript line (gene_id and txpt_id info). Skipping transcript. %s\n", buf);
        while (fgets(buf, 100000, fp)!=NULL) {
          if (strstr(buf, "\ttranscript\t"))
            goto next_transcript;
        }
        return;  /* reached end of file */
      }

      if ((idxc=lookup_chrom(chrom, Chroms, numChroms))<0) {
        fprintf(stderr, "Warning: cannot find chromosome. Skipping transcript. %s\n", buf);
        while (fgets(buf, 100000, fp)!=NULL) {
          if (strstr(buf, "\ttranscript\t"))
            goto next_transcript;
        }
        return;  /* reached end of file */
      }

      if (*numtranscripts>=*maxtranscripts) {
        *maxtranscripts = max(DEFAULT_ALLOC, 2*(*maxtranscripts));
        *transcripts = realloc(*transcripts, (*maxtranscripts)*sizeof(transcript_t));
      }

      txpt = (*transcripts)+(*numtranscripts);

      txpt->gaid = idxc;
      txpt->ori = ori;
      strcpy(txpt->gene_id,sdummy1);
      strcpy(txpt->transcript_id,sdummy2);
      strcpy(txpt->gene_name,sdummy3);
      txpt->fpkm = fpkm;

      txpt->num_exons = txpt->max_exons = txpt->len = 0;
      txpt->from = txpt->to = NULL;
      txpt->printed = 0;

      /* now read the exons */
      while (fgets(buf, 100000, fp)!=NULL) {
        if (strstr(buf, "\ttranscript\t")) {

          txpt->text = ckalloc(txpt->num_exons*(9*2+2)*sizeof(char));
          transcript2text(txpt, txpt->text);
//        printf("%s\n",txpt->text);

          (*numtranscripts)++; 

          goto next_transcript;
        } 

        if (txpt->num_exons>=txpt->max_exons) {
          txpt->max_exons = max(DEFAULT_ALLOC_SMALL, 2*txpt->max_exons);
          txpt->from = realloc(txpt->from, txpt->max_exons*sizeof(int));
          txpt->to = realloc(txpt->to, txpt->max_exons*sizeof(int));
        }

        if (sscanf(buf, "%s\tCufflinks\texon\t%d\t%d\t%d\t%c\t.\t",
                       chrom, &from, &to, &score, &ori)!=5) {
          fprintf(stderr, "Warning: Unrecognized exon line. Skipping transcript. %s\n", buf);

          txpt->num_exons = 0; 
          while (fgets(buf, 100000, fp)!=NULL) {
            if (strstr(buf, "\ttranscript\t"))
              goto next_transcript;
          }
          return;  /* reached the end of file */
        }
        txpt->from[txpt->num_exons] = from;
        txpt->to[txpt->num_exons] = to;
        txpt->len += to-from+1;

        txpt->num_exons++;
      }

      txpt->text = ckalloc(txpt->num_exons*(9*2+2)*sizeof(char));
      transcript2text(txpt, txpt->text);

      (*numtranscripts)++;
   }

   fclose(fp);
}
#endif

/* event_id        event_type      gene_id chrom   event_start     event_end       event_pattern   strand
   1000001 TSS     ENSMUSG00000001138      chr1    36568705        36569998        36569998        +
*/
static int readEventLine(char *buf, event_t *eptr)
{
   char *tok;

   if ((tok = strtok(buf, "\t\n"))==NULL) return 1; 
   eptr->event_id = atoi(tok);

   if ((tok = strtok(NULL, "\t\n"))==NULL) return 1; 
   eptr->type = lookup_event_type(tok);
   if (eptr->type < 0) return 1;
    
   if ((tok = strtok(NULL, "\t\n"))==NULL) return 1; 
   strcpy(eptr->gene_id, tok);

   if ((tok = strtok(NULL, "\t\n"))==NULL) return 1; 
   eptr->gaid = lookup_chrom(tok, Chroms, numChroms);
   if (eptr->gaid < 0) return 1;

   if ((tok = strtok(NULL, "\t\n"))==NULL) return 1; 
   eptr->from = atoi(tok);
   if (eptr->from == 0) return 1;

   if ((tok = strtok(NULL, "\t\n"))==NULL) return 1; 
   eptr->to = atoi(tok);
   if (eptr->to == 0) return 1;

   if ((tok = strtok(NULL, "\t\n"))==NULL) return 1; 
   strcpy(eptr->pattern, tok);

   if ((tok = strtok(NULL, "\t\n"))==NULL) return 1; 
   eptr->ori = *tok;
   if ((eptr->ori!='+') && (eptr->ori!='-') && (eptr->ori!='.')) return 1;

   eptr->fpkm = 0;

   memset(eptr->other,0,100);

   return 0;
}


/****************      Utilities      ****************/
static int lookup_event_type(char *subj)
{
   int i;

   for (i=0; (i<numCodes) && strcmp(Codes[i],subj); i++);
   return (i==numCodes) ? -1 : i;
}

static int lookup_chrom(char *subj, char *chroms[], int numchroms)
{
   int i;

   for (i=0; (i<numchroms) && strcmp(chroms[i],subj); i++);
   return (i==numChroms) ? -1 : i;
}

static int  transcript_cmp(const void *a, const void *b) {
   transcript_t *A = (transcript_t *)a;
   transcript_t *B = (transcript_t *)b;
   int     t;

   if ((t=(A->gaid-B->gaid))!=0) return t;

   return (A->from[0]-B->from[0]);
}

static int overlap(int a, int b, int x, int y)
{
  return ((a<=x && x<=b) || (x<=a && a<=y));
}

static int  searchExonApprox(int f, int t, transcript_t *tptr, int W, int *xi)
{
  int i;

  *xi = -1;

  for (i=0; i<tptr->num_exons && tptr->from[i]<=t; i++)
    if (abs(tptr->from[i]-f)<=W && abs(tptr->to[i]-t)<=W) {
      *xi = i;
      return 0; 
  }
   
  return 1;
}

static int  searchIntronApprox(int f, int t,transcript_t *tptr, int W, int *yi)
{
  int i;

  *yi = -1;

  for (i=0; i<tptr->num_exons-1; i++)
    if (abs(tptr->to[i]-f)<=W && abs(tptr->from[i+1]-t)<=W) {
      *yi = i+1;
      return 0;
    }

  return 1;
}

static void loadJuncs(const char *filename, junc_t **juncs, int *numjuncs, int *maxjuncs)
{
   FILE         *fp = ckopen(filename, "r");
   char          buf[100000], chrom[20], ori;
   int           from, to, num, idxc;
   junc_t       *jptr;

   *numjuncs = *maxjuncs = 0;

   while (fgets(buf, 100000, fp)!=NULL) {

      if (buf[0] == '#') continue;

      /* chr1 14829 14970 150 - */

      if (sscanf(buf, "%s\t%d\t%d\t%d\t%c",
                       chrom, &from, &to, &num, &ori) != 5) {
        fprintf(stderr, "Warning: Unrecognized junction line. Skipping junction. %s\n", buf);
        continue;
      }

      if ((idxc=lookup_chrom(chrom, Chroms, numChroms))<0) {
        fprintf(stderr, "Warning: cannot find chromosome. Skipping junction. %s\n", buf);
        continue;
      }

      if (*numjuncs>=*maxjuncs) {
        *maxjuncs = max(DEFAULT_ALLOC, 2*(*maxjuncs));
        *juncs = realloc(*juncs,(*maxjuncs)*sizeof(junc_t));
      }

      jptr = (*juncs)+(*numjuncs);

      jptr->from = from;
      jptr->to = to;
      jptr->num = num;
      jptr->gaid = idxc;
      jptr->ori = ori;

      (*numjuncs)++;
   }

   fclose(fp);
}

static int  junc_cmp(const void *a, const void *b) {
   junc_t *A = (junc_t *)a;
   junc_t *B = (junc_t *)b;
   int     t;

   if ((t=(A->gaid-B->gaid))!=0) return t;

   return (A->from-B->from);
}

static int searchJuncApprox(int il, int ir, int W, junc_t *juncs, int numjuncs)
{
   int num = 0, idx;

   idx = junc_lookup2(il-W,juncs,0,numjuncs-1);
   while (idx>0 && idx<numjuncs && abs(juncs[idx].from-il)<=W) {
     if (abs(juncs[idx].to-ir)<=W) num += juncs[idx].num;
     idx++;
   }

   return num; 
}
   
/* return the first entry larger or equal to the value x */
static int junc_lookup2(int x, junc_t *juncs, int l, int u)
{
   int m = (int)((l+u)/2);

   if (x <= juncs[l].from) return l;

   if (x > juncs[u].from) return -1;

   if ((u-l<=1) && (x<=juncs[u].from)) return u;
     
#ifdef DEBUG
   printf("%d %d(%d) %d(%d) %d(%d)\n", x, l, juncs[l].from, m, juncs[m].from, u, juncs[u].from);
#endif

   assert(m+1<=u);
   assert(m-1>=l);

   if (juncs[m].from < x) return junc_lookup2(x,juncs,m+1,u);
  
   if (juncs[m].from>=x) {
     if (juncs[m-1].from<x) return m;
     else return junc_lookup2(x,juncs,l,m-1);
   }
}

