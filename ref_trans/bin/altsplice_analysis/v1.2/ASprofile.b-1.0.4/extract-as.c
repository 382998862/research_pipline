/*
   extract-as:

   Software for discovery of within-gene splicing variation from
   transcripts produced from RNA-seq data, and comparison to
   a reference annotation:

   Types of variation include skipped exons or exon cassettes,   
   and alternative transcription start and termination:
       - TSS alternative transcription start (5' terminal exon)
       - TTS alternative transcription termination (3' terminal end)
       - SKIP(_ON,_OFF) single exon skip, exact intron boundaries
       - XSKIP(_ON,_OFF) single exon skip, approximate/variable intron boundaries
       - MSKIP(_ON,_OFF) multiple exon skip, exact splice boundaries
       - XMSKIP(_ON,_OFF) multiple exon skip, approximate/variable intron boundaries
       - IR(_ON,_OFF) single intron retained, exact splice boundaries
       - XIR(_ON,_OFF) single intron retained, approximate/variable exon boundaries 
       - MIR(_ON,_OFF) multiple introns retained, exact exon boundaries
       - XMIR(_ON,_OFF) multiple introns retained, approximate/variable exon boundaries 
       - AE alternative exon end(s), flanking intron matches exactly at the distal end 
       - XAE alternative exon end(s), flanking intron matches approximately at the distal end

   Copyright (C) 2011-2013  Liliana Florea (GPL)

   Compile with: gcc extract-as.c lib.c -o extract_as -lm
*/

#include "libc.h"
#include "lib.h"

#define DEFAULT_ALLOC        10000
#define DEFAULT_ALLOC_SMALL    100

#define V      100
#define W       10

static char Usage[] = "%s txptgtf genome_hdrs [-r tmap goldref]\n";

static char **Chroms;
static int    numChroms = 0;
static int    maxChroms = 0;

typedef struct {
   char transcript_id[50];
   char gene_id[50];
   char gene_name[50];
   int  gaid;
   char ori;
   int  num_exons;
   int  max_exons;
   int *from, *to;
   double fpkm;
   int  printed;
   int  len;
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

typedef struct event {
  int  event_id;
  int  type;
  int  gaid, from, to;
  char pattern[10000], ori;
} event_t;

typedef struct tmap {
  char txpt_id[100], gene_id[100];
  char goldref_id[100];
  char class_code;
} tmap_t;


static char      errmsg[100];

static int       eventCounter = 1000000;
static event_t  *Q=NULL;
static int       nQ=0, maxQ=0;

static void  extract_as(transcript_t *,int,int,tmap_t *,int,transcript_t *,int);
static void  extract_as_tss(transcript_t *,int,int,tmap_t *,int,transcript_t *,int);
static void  extract_as_tts(transcript_t *,int,int,tmap_t *,int,transcript_t *,int);
static void  extract_as_skipped(transcript_t *,int,int,tmap_t *,int,transcript_t *,int);
static void  extract_as_skippedX(transcript_t *,int,int,tmap_t *,int,transcript_t *,int);
static void  extract_as_ir(transcript_t *,int,int,tmap_t *,int,transcript_t *,int);
static void  extract_as_ae(transcript_t *,int,int,tmap_t *,int,transcript_t *,int);
static void  load_chroms(char *gafile, char ***chroms, int *numchroms, int *maxchroms);
#if 0
static void  loadTranscripts(const char *,transcript_t **,int *,int *);
#endif
static int   transcript_cmp(const void *a,const void *b);
static int   goldref_cmp(const void *a,const void *b);
static void  loadTmap(const char *,tmap_t **,int *,int *);
static void  loadGTF(const char *,transcript_t **,int *,int *);
static void  init_event(void);
static int   event_lookup_and_add(int gaid,int from,int to,char ori,int type,char *pattern);
static int   lookup_chrom(char *subj, char *chroms[], int numchroms);
static int   lookup_tmap(char *transcript_id, tmap_t *tmaps, int numtmaps);
static int   lookup_tmap2(char *transcript_id, tmap_t *tmaps, int l, int u);
static int   lookup_goldref(char *goldref_id, transcript_t *goldrefs, int numgoldrefs);
static int   lookup_goldref2(char *goldref_id, transcript_t *goldrefs, int l, int u);
static int   alternative_tss(transcript_t *t1, transcript_t *t2, int type);
static void  exon_match(int from, int to, transcript_t *g, char *sign, int *xo);
static int   overlap(int a, int b, int x, int y);
inline void  print_goldref_context(char *,int,int,tmap_t *,int,transcript_t *,int);

inline void  print_goldref_context(char *name, int from, int to, tmap_t *tmaps, int numtmaps, transcript_t *goldrefs, int numgoldrefs)
{
   int idxm, idxg;

   if (goldrefs == NULL) printf("\n");
   else {
     idxm = lookup_tmap(name, tmaps, numtmaps);
     if (idxm<=0)
       printf("\t-\t-\t-\n");
     else {
       char sign; int x;

       idxg = lookup_goldref(tmaps[idxm].goldref_id, goldrefs, numgoldrefs);
       assert(idxg>=0);

       exon_match(from, to, goldrefs+idxg, &sign, &x);
       printf("\t%s\t%c\t%c%d\n", tmaps[idxm].goldref_id, tmaps[idxm].class_code, sign, x);
     }
   }
}

static int  transcript_cmp(const void *a, const void *b) {
   transcript_t *A = (transcript_t *)a;
   transcript_t *B = (transcript_t *)b;
   int     t;

   if ((t=(A->gaid-B->gaid))!=0) return t;

   return strcmp(A->gene_id,B->gene_id);
}

static int  goldref_cmp(const void *a, const void *b) {
   transcript_t *A = (transcript_t *)a;
   transcript_t *B = (transcript_t *)b;
   int     t;

   if ((t=(A->gaid-B->gaid))!=0) return t;

   return strcmp(A->transcript_id,B->transcript_id);
}


static int  tmap_cmp(const void *a, const void *b) {
   tmap_t *A = (tmap_t *)a;
   tmap_t *B = (tmap_t *)b;

   return strcmp(A->txpt_id,B->txpt_id);
}


void main(int argc, char *argv[])
{
    transcript_t **Transcripts;
    int           *transcriptCounts;
    int            numTranscripts, maxTranscripts;
    transcript_t **Goldref = NULL;
    int           *GoldrefCounts = NULL;
    int            numGoldref, maxGoldref;
    tmap_t        *Tmap;
    int            numTmap, maxTmap;
    int            start, end, i;
  
    if (argc!=3 && argc!=6) fatalf(Usage, argv[0]);

    if (argc==6 && strcmp(argv[3], "-r")) fatalf(Usage, argv[0]);
 
    load_chroms(argv[2], &Chroms, &numChroms, &maxChroms);

    numGoldref = maxGoldref = 0; Goldref = NULL;
    numTmap = maxTmap = 0; Tmap = NULL;

    if (argc==6) {

      loadTmap(argv[4], &Tmap, &numTmap, &maxTmap);

      qsort(Tmap, numTmap, sizeof(tmap_t), tmap_cmp); /* alphabetically by cufflinks gene name */

      Goldref = (transcript_t **)malloc(numChroms*sizeof(transcript_t));
      memset(Goldref, 0, numChroms*sizeof(transcript_t *));

      loadGTF(argv[5], &(Goldref[0]), &numGoldref, &maxGoldref);

      qsort(Goldref[0], numGoldref, sizeof(transcript_t), goldref_cmp);
      GoldrefCounts = (int *)malloc(numChroms*sizeof(int));
      memset(GoldrefCounts, 0, numChroms*sizeof(int));
   
      for (i=0; i<numGoldref; i++) {
        int gaid = Goldref[0][i].gaid;
        if (Goldref[gaid]==NULL) Goldref[gaid] = Goldref[0]+i;
        GoldrefCounts[gaid]++;
      }
    }

    numTranscripts = maxTranscripts = 0;

    Transcripts = (transcript_t **)malloc(numChroms*sizeof(transcript_t));
    memset(Transcripts, 0, numChroms*sizeof(transcript_t *));

    loadGTF(argv[1], &(Transcripts[0]), &numTranscripts, &maxTranscripts);
    qsort(Transcripts[0], numTranscripts, sizeof(transcript_t), transcript_cmp); 
    transcriptCounts = (int *)malloc(numChroms*sizeof(int));
    memset(transcriptCounts, 0, numChroms*sizeof(int));

    for (i=0; i<numTranscripts; i++) {
      int gaid = Transcripts[0][i].gaid;
      if (Transcripts[gaid]==NULL) Transcripts[gaid] = Transcripts[0]+i; 
      transcriptCounts[gaid]++;
    }

    /* print header */
    printf("#event_id\tevent_type\tgene_id\tchrom\tevent_start\tevent_end\tevent_pattern\tstrand");
    printf("\ttxpt_id\ttxpt_numexons\tthis_exon\tnum_alt_exons\texon_froms\texon_tos\ttxpt_fpkm\ttxpt_len");

    if (Goldref==NULL) 
      printf("\n");
    else 
      printf("\tref_id\tclass_code\tref_exon\n");

    for (i=start=end=0; i<numTranscripts; i++) {
      if (strcmp(Transcripts[0][i].gene_id,Transcripts[0][start].gene_id)) {
        int gaid = Transcripts[0][start].gaid;
        extract_as(Transcripts[0], start, end, Tmap, numTmap,
                   (Goldref!=NULL) ? Goldref[gaid] : NULL,
                   (Goldref!=NULL) ? GoldrefCounts[gaid] : 0);
        start = i;
      }
      end = i;
    }
    if (numTranscripts) {
      int gaid = Transcripts[0][numTranscripts-1].gaid;
      extract_as(Transcripts[0], start, end, Tmap, numTmap,
                 (Goldref!=NULL) ? Goldref[gaid] : NULL,
                 (Goldref!=NULL) ? GoldrefCounts[gaid] : 0);
    }

    free(Transcripts[0]);
    free(Transcripts);

    free(Tmap);

    if (Goldref != NULL) {
      free(Goldref[0]);
      free(Goldref);
    }

    exit(0);
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
          strncpy(sdummy3, tok+1, strlen(tok)-1);
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
      (*numtranscripts)++;
   }
   
   fclose(fp);
}
#endif

static void loadTmap(const char *filename, tmap_t **tmaps, int *numtmaps, int *maxtmaps)
{
   FILE    *fp = ckopen(filename, "r");
   char     buf[100000], *tok;
   tmap_t  *t;

   *numtmaps = *maxtmaps = 0;

   fgets(buf, 100000, fp);   /* skip over header line */

   while (fgets(buf, 100000, fp)!=NULL) {
     if (buf[0] == '-') continue;

     if (*numtmaps>=*maxtmaps) {
       *maxtmaps = max(DEFAULT_ALLOC_SMALL, 2*(*maxtmaps));
       *tmaps = realloc(*tmaps, (*maxtmaps)*sizeof(tmap_t));
     }
     t = (*tmaps)+(*numtmaps);

     tok = strtok(buf, "\t\n");  /* read and skip over goldref gene_id */
     tok = strtok(NULL, "\t\n");
     strcpy(t->goldref_id, tok);
     tok = strtok(NULL, "\t\n");
     t->class_code = tok[0];
     tok = strtok(NULL, "\t\n");
     strcpy(t->gene_id, tok);
     tok = strtok(NULL, "\t\n");
     strcpy(t->txpt_id, tok);

     (*numtmaps)++;
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

     if ((txpt->num_exons>1) && (txpt->ori=='-') && (txpt->from[0]>txpt->from[1])) {
       int u, l, aux;

       l = 0; u = txpt->num_exons-1;
       while (l<u) {
         aux = txpt->from[l]; txpt->from[l] = txpt->from[u]; txpt->from[u] = aux;
         aux = txpt->to[l]; txpt->to[l] = txpt->to[u]; txpt->to[u] = aux;

         l++; u--;
       }
     }
   }
}


/* TSS/TTS: it will only print if there is an alternative splicing event there;
        single TSS(TTS) genes are not reported. May print multiple (incomplete) starts
        for a TSS(TTS), but with the same event_id. These can be filtered later.
*/
static void extract_as(transcript_t *transcripts, int start, int end, tmap_t *tmaps, int numtmaps, transcript_t *goldrefs, int numgoldrefs)
{

   int i;

   for (i=start; i<=end; i++)
     if (transcripts[i].ori != '.') break;

   if (i>end) return;

   /* find alternative TSSs */
   
   init_event();

   extract_as_tss(transcripts,start,end,tmaps,numtmaps,goldrefs,numgoldrefs);

   extract_as_tts(transcripts,start,end,tmaps,numtmaps,goldrefs,numgoldrefs);

   extract_as_skippedX(transcripts,start,end,tmaps,numtmaps,goldrefs,numgoldrefs);

   extract_as_ir(transcripts,start,end,tmaps,numtmaps,goldrefs,numgoldrefs);

   extract_as_ae(transcripts,start,end,tmaps,numtmaps,goldrefs,numgoldrefs);

   return;
}

static void extract_as_tss(transcript_t *transcripts, int start, int end, tmap_t *tmaps, int numtmaps,transcript_t *goldrefs, int numgoldrefs)
{
   int           i, i0, j, k, idx, hasAS;
   char          pattern[10000];
   transcript_t *t1, *t2;

   hasAS = 0;

   for (i=start; i<end; i++) {
     t1 = transcripts+i;
     if (t1->num_exons == 1) continue;

     for (j=i+1; j<=end; j++) {
       t2 = transcripts+j;
       if (t2->num_exons == 1) continue;

       if (alternative_tss(t1,t2,TSS_TYPE)) {
         /* have one internal exon edge in common and different starts */
         if (t1->printed == 0) {
           int from, to;

           if (t1->ori == '+') {
             from = t1->from[0]; to = t1->to[0];
           } else {
             from = t1->from[t1->num_exons-1]; to = t1->to[t1->num_exons-1];
           }

           sprintf(pattern, "%d", (t1->ori=='+') ? to : from);     /* match only the 'fixed' (known) end, depending on ori */
           idx = event_lookup_and_add(t1->gaid,from,to,t1->ori,TSS_TYPE,pattern);

           /* this_exon = 1 */
           printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t1\t1\t",
                   Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid], from, to,
                   Q[idx].pattern, t1->ori, t1->transcript_id, t1->num_exons);
           printf("%d",t1->from[0]);
           for (k=1; k<t1->num_exons; k++) printf(",%d", t1->from[k]);
           printf("\t%d",t1->to[0]);
           for (k=1; k<t1->num_exons; k++) printf(",%d", t1->to[k]);
           printf("\t%f\t%d", t1->fpkm, t1->len);

           print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

           t1->printed = 1;

           hasAS = 1;
         }
         if (t2->printed == 0) {
           int from, to;

           if (t2->ori == '+') {
             from = t2->from[0]; to = t2->to[0];
           } else {
             from = t2->from[t2->num_exons-1]; to = t2->to[t2->num_exons-1];
           }

           sprintf(pattern, "%d", (t2->ori=='+') ? to : from);     /* match only the 'fixed' (known) end */
           idx = event_lookup_and_add(t2->gaid,from,to,t2->ori,TSS_TYPE,pattern);

           /* this_exon = 1 */
           printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t1\t1\t",
                   Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid], from, to,
                   Q[idx].pattern, t2->ori, t2->transcript_id, t2->num_exons);
           printf("%d",t2->from[0]);
           for (k=1; k<t2->num_exons; k++) printf(",%d", t2->from[k]);
           printf("\t%d",t2->to[0]);
           for (k=1; k<t2->num_exons; k++) printf(",%d", t2->to[k]);
           printf("\t%f\t%d", t2->fpkm, t2->len);

           print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

           t2->printed = 1;
  
           hasAS = 1;
         }
       }  /* alternative_tss */
     }  /* for j */
   }  /* for i */

#if 1
   if (hasAS == 0) {
     int from, to;

     for (i=start, from=to=i0=-1; i<=end; i++) {
       t1 = transcripts+i;
       if (t1->num_exons == 1) continue;

       if (t1->ori == '+') {
         if ((to<0) || (to>t1->to[0])) {
           from = t1->from[0]; to = t1->to[0]; i0 = i;
         }
       } else {
         if ((from<0) || (from<t1->to[t1->num_exons-1])) {
           from = t1->from[t1->num_exons-1]; to = t1->to[t1->num_exons-1]; i0 = i;
         }
       }
     }

     if (i0<0) return;   /* single exon transcripts only; may not know orientation, so return */

     t1 = transcripts+i0;
     sprintf(pattern, "%d", (t1->ori=='+') ? to : from);     /* match only the 'fixed' (known) end, depending on ori */
     idx = event_lookup_and_add(t1->gaid,from,to,t1->ori,TSS_TYPE,pattern);

     /* this_exon = 1 */
     printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t1\t1\t",
            Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid], from, to,
            Q[idx].pattern, t1->ori, t1->transcript_id, t1->num_exons);
     printf("%d",t1->from[0]);
     for (k=1; k<t1->num_exons; k++) printf(",%d", t1->from[k]);
     printf("\t%d",t1->to[0]);
     for (k=1; k<t1->num_exons; k++) printf(",%d", t1->to[k]);
     printf("\t%f\t%d", t1->fpkm, t1->len);

     print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);
   }
#endif
}

static void extract_as_tts(transcript_t *transcripts, int start, int end, tmap_t *tmaps, int numtmaps, transcript_t *goldrefs, int numgoldrefs)
{
   int           i, i0, j, k, idx, hasAS;
   char          pattern[10000];
   transcript_t *t1, *t2; 

   hasAS = 0;

   for (i=start; i<end; i++) {
     t1 = transcripts+i;
     if (t1->num_exons == 1) continue;

     for (j=i+1; j<=end; j++) {
       t2 = transcripts+j;
       if (t2->num_exons == 1) continue;

       if (alternative_tss(t1,t2,TTS_TYPE)) {
         /* have one internal exon edge in common and different starts */
         if (t1->printed == 0) {
           int from, to;

           if (t1->ori == '+') {
             from = t1->from[t1->num_exons-1]; to = t1->to[t1->num_exons-1];
           } else {
             from = t1->from[0]; to = t1->to[0];
           }

           sprintf(pattern, "%d", (t1->ori=='+') ? from : to);     /* match only the 'fixed' (known) end, depending on ori */
           idx = event_lookup_and_add(t1->gaid,from,to,t1->ori,TTS_TYPE,pattern);

           /* this_exon = 1 */
           printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t1\t",
                   Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid], from, to,
                   Q[idx].pattern, t1->ori, t1->transcript_id, t1->num_exons, t1->num_exons);
           printf("%d",t1->from[0]);
           for (k=1; k<t1->num_exons; k++) printf(",%d", t1->from[k]);
           printf("\t%d",t1->to[0]);
           for (k=1; k<t1->num_exons; k++) printf(",%d", t1->to[k]);
           printf("\t%f\t%d", t1->fpkm, t1->len);

           print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

           t1->printed = 1;

           hasAS = 1;
         }
         if (t2->printed == 0) {
           int from, to;

           if (t2->ori == '+') {
             from = t2->from[t2->num_exons-1]; to = t2->to[t2->num_exons-1];
           } else {
             from = t2->from[0]; to = t2->to[0];
           }

           sprintf(pattern, "%d", (t2->ori=='+') ? from : to);     /* match only the 'fixed' (known) end */
           idx = event_lookup_and_add(t2->gaid,from,to,t2->ori,TTS_TYPE,pattern);

           /* this_exon = 1 */
           printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t1\t",
                   Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid], from, to,
                   Q[idx].pattern, t2->ori, t2->transcript_id, t2->num_exons, t2->num_exons);
           printf("%d",t2->from[0]);
           for (k=1; k<t2->num_exons; k++) printf(",%d", t2->from[k]);
           printf("\t%d",t2->to[0]);
           for (k=1; k<t2->num_exons; k++) printf(",%d", t2->to[k]);
           printf("\t%f\t%d", t2->fpkm, t2->len);

           print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

           t2->printed = 1;

           hasAS = 1;
         }
       }  /* alternative_tts */
     }  /* for j */
   }  /* for i */

#if 1
   if (hasAS == 0) {
     int from, to;

     for (i=start, from=to=i0=-1; i<=end; i++) {
       t1 = transcripts+i;
       if (t1->num_exons == 1) continue;

       if (t1->ori == '+') {
         if ((from<0) || (from<t1->from[t1->num_exons-1])) {
           from = t1->from[t1->num_exons-1]; to = t1->to[t1->num_exons-1]; i0 = i;
         }
       } else {
         if ((to<0) || (to>t1->to[0])) {
           from = t1->from[0]; to = t1->to[0]; i0 = i;
         }
       }
     }

     if (i0<0) return;   /* single exon transcripts only; may not know orientation, so return */
 
     t1 = transcripts+i0;

     sprintf(pattern, "%d", (t1->ori=='+') ? from : to);     /* match only the 'fixed' (known) end, depending on ori */
     idx = event_lookup_and_add(t1->gaid,from,to,t1->ori,TTS_TYPE,pattern);

     /* this_exon = 1 */
     printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t1\t",
            Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid], from, to,
            Q[idx].pattern, t1->ori, t1->transcript_id, t1->num_exons, t1->num_exons);
     printf("%d",t1->from[0]);
     for (k=1; k<t1->num_exons; k++) printf(",%d", t1->from[k]);
     printf("\t%d",t1->to[0]);
     for (k=1; k<t1->num_exons; k++) printf(",%d", t1->to[k]);
     printf("\t%f\t%d", t1->fpkm, t1->len);

     print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);
   }
#endif
}

static void extract_as_skipped(transcript_t *transcripts, int start, int end, tmap_t *tmaps, int numtmaps, transcript_t *goldrefs, int numgoldrefs)
{
   int           i, j, k, l, u, v, ui, vi, x, idx;
   char          pattern[10000], *s;
   transcript_t *t1, *t2;
   
   for (i=start; i<end; i++) {
     t1 = transcripts+i;
     if (t1->num_exons<2) continue;

     for (j=i+1; j<=end; j++) {
       t2 = transcripts+j;
       if (t2->num_exons<2) continue;

       if (t1->num_exons+t2->num_exons<5) continue;

       k = l = 0;
       while ((k<t1->num_exons-1) && (l<t2->num_exons-1)) {
         if (t1->to[k]<t2->to[l]) k++;
         else if (t2->to[l]<t1->to[k]) l++;
         else { /* equal; check for alternative splicing */

           if (t1->from[k+1]<t2->from[l+1]) {
             for (u=1; (k+u<t1->num_exons) && (t1->from[k+u]<t2->from[l+1]); u++);

             if ((u>1) && (k<t1->num_exons-u) && (t1->from[k+u]==t2->from[l+1])) {

               /* t1 'on'; t2 'off' */

               int from, to, type;

               if (t1->ori=='+') {
                 from = t1->from[k+1]; to = t1->to[k+1];
               } else {
                 from = t1->from[k+u-1]; to = t1->to[k+u-1];
               }

               type = (u>2) ? MSKIP_ON_TYPE : SKIP_ON_TYPE;

               memset(pattern, 0, 10000); s = pattern;
               sprintf(s, "%d,%d-%d", t1->to[k], t1->from[k+1], t1->to[k+1]);     /* match whole exon(s) pattern */
               while (*s) s++;
               for (ui=2; ui<u; ui++) {
                 sprintf(s, ",%d-%d", t1->from[k+ui], t1->to[k+ui]);
                 while (*s) s++;
               }
               sprintf(s, ",%d", t1->from[k+u]);
               idx = event_lookup_and_add(t1->gaid,t1->from[k+1],t1->to[k+u-1],t1->ori,type,pattern);

               printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                      Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                      t1->from[k+1], t1->to[k+u-1], Q[idx].pattern, t1->ori,
                      t1->transcript_id, t1->num_exons,
                      ((t1->ori=='+') ? (k+2) : (t1->num_exons-(k+u-1))), u-1);
               printf("%d",t1->from[0]);
               for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
               printf("\t%d",t1->to[0]);
               for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
               printf("\t%f\t%d", t1->fpkm, t1->len);

               print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);


               /* now the 'off' form */
               type = (u>2) ? MSKIP_OFF_TYPE : SKIP_OFF_TYPE;

               sprintf(pattern, "%d,%d", t2->to[l], t2->from[l+1]);  /* surrounding exon ends */
               idx = event_lookup_and_add(t1->gaid,t1->from[k+1],t1->to[k+u-1],t1->ori,type,pattern);

               printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t>%d\t%d\t",
                      Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                      t1->from[k+1], t1->to[k+u-1], Q[idx].pattern, t2->ori,
                      t2->transcript_id, t2->num_exons,
                      ((t2->ori == '+') ? (l+1) : (t2->num_exons-(l+1))), (u-1));
               printf("%d",t2->from[0]);
               for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
               printf("\t%d",t2->to[0]);
               for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
               printf("\t%f\t%d", t2->fpkm, t2->len);

               print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);


               k += u; l += 1;
             } else {
               k++; l++;
             }
           } else if (t2->from[l+1]<t1->from[k+1]) {
             for (v=1; (l+v<t2->num_exons) && (t2->from[l+v]<t1->from[k+1]); v++);

             if ((v>1) && (l<t2->num_exons-v) && (t2->from[l+v]==t1->from[k+1])) {

               /* t1 'off'; t2 'on' */

               int from, to, type;

               if (t2->ori=='+') {
                 from = t2->from[l+1]; to = t2->to[l+1];
               } else {
                 from = t2->from[l+v-1]; to = t2->to[l+v-1];
               }

               type = (v>2) ? MSKIP_ON_TYPE : SKIP_ON_TYPE;

               memset(pattern, 0, 10000); s = pattern;
               sprintf(s, "%d,%d-%d", t2->to[l], t2->from[l+1], t2->to[l+1]);     /* match whole exon(s) pattern */
               while (*s) s++;
               for (vi=2; vi<v; vi++) {
                 sprintf(s, ",%d-%d", t2->from[l+vi], t2->to[l+vi]);
                 while (*s) s++;
               }
               sprintf(s, ",%d", t2->from[l+v]);
               idx = event_lookup_and_add(t2->gaid,t2->from[l+1],t2->to[l+v-1],t2->ori,type,pattern);

               printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                      Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                      t2->from[l+1], t2->to[l+v-1], Q[idx].pattern, t2->ori,
                      t2->transcript_id, t2->num_exons,
                      ((t2->ori=='+') ? (l+2) : (t2->num_exons-(l+v-1))), (v-1));
               printf("%d",t2->from[0]);
               for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
               printf("\t%d",t2->to[0]);
               for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
               printf("\t%f\t%d", t2->fpkm, t2->len);

               print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);


               /* now the 'off' form */
               type = (v>2) ? MSKIP_OFF_TYPE : SKIP_OFF_TYPE;

               sprintf(pattern, "%d,%d", t1->to[k], t1->from[k+1]);  /* surrounding exon ends */
               idx = event_lookup_and_add(t1->gaid,t1->to[k],t1->from[k+1],t1->ori,type,pattern);

               printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t>%d\t%d\t",
                      Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                      t2->from[l+1], t2->to[l+v-1], Q[idx].pattern, t1->ori,
                      t1->transcript_id, t1->num_exons,
                      ((t1->ori == '+') ? (k+1) : (t1->num_exons-(k+1))), (v-1));
               printf("%d",t1->from[0]);
               for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
               printf("\t%d",t1->to[0]);
               for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
               printf("\t%f\t%d", t1->fpkm, t1->len);

               print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);


               l += v; k += 1;
             } else {
               k++; l++;
             }
           } else {
             /* no alternative skipped exon here */
             k++; l++;
           }
         } /* end equal to's */
       } /* end criss-cross */
     } /* for j */
   } /* for i */

}


static void extract_as_skippedX(transcript_t *transcripts, int start, int end, tmap_t *tmaps, int numtmaps, transcript_t *goldrefs, int numgoldrefs)
{
   int           i, j, k, l, u, v, ui, vi, x, idx;
   char          pattern[10000], *s;
   transcript_t *t1, *t2;

   for (i=start; i<end; i++) {
     t1 = transcripts+i;
     if (t1->num_exons<2) continue;

     for (j=i+1; j<=end; j++) {
       t2 = transcripts+j;
       if (t2->num_exons<2) continue;

       if (t1->num_exons+t2->num_exons<5) continue;

       k = l = 0;
       while ((k<t1->num_exons-1) && (l<t2->num_exons-1)) {
         if (t1->to[k]<t2->from[l]) k++;
         else if (t2->to[l]<t1->from[k]) l++;
         else { /* overlap; check for alternative splicing */

           if (t1->from[k+1]<t2->from[l+1]) {
             for (u=1; (k+u<t1->num_exons) && (t2->to[l]<t1->from[k+u]) && (t1->to[k+u]<t2->from[l+1]); u++);

             if ((u>1) && (k+u<t1->num_exons) && overlap(t1->from[k+u],t1->to[k+u],t2->from[l+1],t2->to[l+1])) {

               /* t1 'on'; t2 'off' */

               int from, to, type;

               if (t1->ori=='+') {
                 from = t1->from[k+1]; to = t1->to[k+1];
               } else {
                 from = t1->from[k+u-1]; to = t1->to[k+u-1];
               }

               if ((t1->to[k]==t2->to[l]) && (t1->from[k+u]==t2->from[l+1]))
                 type = (u>2) ? MSKIP_ON_TYPE : SKIP_ON_TYPE;
               else
                 type = (u>2) ? XMSKIP_ON_TYPE : XSKIP_ON_TYPE;

               memset(pattern, 0, 10000); s = pattern;
               sprintf(s, "%d,%d-%d", t1->to[k], t1->from[k+1], t1->to[k+1]);     /* match whole exon(s) pattern */
               while (*s) s++;
               for (ui=2; ui<u; ui++) {
                 sprintf(s, ",%d-%d", t1->from[k+ui], t1->to[k+ui]);
                 while (*s) s++;
               }
               sprintf(s, ",%d", t1->from[k+u]);
               idx = event_lookup_and_add(t1->gaid,t1->from[k+1],t1->to[k+u-1],t1->ori,type,pattern);

               printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                      Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                      t1->from[k+1], t1->to[k+u-1], Q[idx].pattern, t1->ori,
                      t1->transcript_id, t1->num_exons,
                      ((t1->ori=='+') ? (k+2) : (t1->num_exons-(k+u-1))), u-1);
               printf("%d",t1->from[0]);
               for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
               printf("\t%d",t1->to[0]);
               for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
               printf("\t%f\t%d", t1->fpkm, t1->len);

               print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

               /* now the 'off' form */
               if ((t1->to[k]==t2->to[l]) && (t1->from[k+u]==t2->from[l+1]))
                 type = (u>2) ? MSKIP_OFF_TYPE : SKIP_OFF_TYPE;
               else
                 type = (u>2) ? XMSKIP_OFF_TYPE : XSKIP_OFF_TYPE;

               sprintf(pattern, "%d,%d", t2->to[l], t2->from[l+1]);  /* surrounding exon ends */
               idx = event_lookup_and_add(t1->gaid,t1->from[k+1],t1->to[k+u-1],t1->ori,type,pattern);

               printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t>%d\t%d\t",
                      Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                      t1->from[k+1], t1->to[k+u-1], Q[idx].pattern, t2->ori,
                      t2->transcript_id, t2->num_exons,
                      ((t2->ori == '+') ? (l+1) : (t2->num_exons-(l+1))), (u-1));
               printf("%d",t2->from[0]);
               for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
               printf("\t%d",t2->to[0]);
               for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
               printf("\t%f\t%d", t2->fpkm, t2->len);

               print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);


               k += u; l += 1;
             } else if (k+u == t1->num_exons) {
               k += u; l += 1;
             } else if (u==1) {
               if (t1->to[k+1]>=t2->from[l+1]) l++; /* overlap(k+1,l+1) */
               k++;
             } else {
               assert(!overlap(t1->from[k+u],t1->to[k+u],t2->from[l+1],t2->to[l+1]));
               k += u; l += 1;
             }
           } else if (t2->from[l+1]<t1->from[k+1]) {

             for (v=1; (l+v<t2->num_exons) && (t1->to[k]<t2->from[l+v]) && (t2->to[l+v]<t1->from[k+1]); v++);

             if ((v>1) && (l+v<t2->num_exons) && overlap(t2->from[l+v],t2->to[l+v],t1->from[k+1],t1->to[k+1])) {

               /* t1 'off'; t2 'on' */

               int from, to, type;

               if (t2->ori=='+') {
                 from = t2->from[l+1]; to = t2->to[l+1];
               } else {
                 from = t2->from[l+v-1]; to = t2->to[l+v-1];
               }

               if ((t2->to[l]==t1->to[k]) && (t1->from[k+1]==t2->from[l+v]))
                 /* classic exon skipping */
                 type = (v>2) ? MSKIP_ON_TYPE : SKIP_ON_TYPE;
               else
                 type = (v>2) ? XMSKIP_ON_TYPE : XSKIP_ON_TYPE;

               memset(pattern, 0, 10000); s = pattern;
               sprintf(pattern, "%d,%d-%d", t2->to[l], t2->from[l+1], t2->to[l+1]);     /* match whole exon(s) pattern */
               while (*s) s++;
               for (vi=2; vi<v; vi++) {
                 sprintf(s, ",%d-%d", t2->from[l+vi], t2->to[l+vi]);
                 while (*s) s++;
               }
               sprintf(s, ",%d", t2->from[l+v]);
               idx = event_lookup_and_add(t2->gaid,t2->from[l+1],t2->to[l+v-1],t2->ori,type,pattern);

               printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                      Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                      t2->from[l+1], t2->to[l+v-1], Q[idx].pattern, t2->ori,
                      t2->transcript_id, t2->num_exons,
                      ((t2->ori=='+') ? (l+2) : (t2->num_exons-(l+v-1))), (v-1));
               printf("%d",t2->from[0]);
               for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
               printf("\t%d",t2->to[0]);
               for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
               printf("\t%f\t%d", t2->fpkm, t2->len);

               print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);


               /* now the 'off' form */
              if ((t2->to[l]==t1->to[k]) && (t1->from[k+1]==t2->from[l+v]))
                 /* classic exon skipping */
                 type = (v>2) ? MSKIP_OFF_TYPE : SKIP_OFF_TYPE;
               else
                 type = (v>2) ? XMSKIP_OFF_TYPE : XSKIP_OFF_TYPE;

               sprintf(pattern, "%d,%d", t1->to[k], t1->from[k+1]);  /* surrounding exon ends */
               idx = event_lookup_and_add(t1->gaid,t1->to[k],t1->from[k+1],t1->ori,type,pattern);

               printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t>%d\t%d\t",
                      Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                      t2->from[l+1], t2->to[l+v-1], Q[idx].pattern, t1->ori,
                      t1->transcript_id, t1->num_exons,
                      ((t1->ori == '+') ? (k+1) : (t1->num_exons-(k+1))), (v-1));
               printf("%d",t1->from[0]);
               for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
               printf("\t%d",t1->to[0]);
               for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
               printf("\t%f\t%d", t1->fpkm, t1->len);

               print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);


               l += v; k += 1;
             } else if (l+v == t2->num_exons) {
               l += v; k += 1;
             } else if (v==1) {
               if (t2->to[l+1]>=t1->from[k+1]) k++; /* overlap(k+1,l+1) */
               l += 1;
             } else {
               assert(!overlap(t2->from[l+v],t2->to[l+v],t1->from[k+1],t1->to[k+1]));
               l += v; k += 1;
             }
           } else {
             /* t1->from[k+1]==t2->from[l+1] */
             k++; l++;
           }
         } /* end equal to's */ 
       } /* end criss-cross */
     } /* for j */ 
   } /* for i */

}

static void extract_as_ir(transcript_t *transcripts, int start, int end, tmap_t *tmaps, int numtmaps, transcript_t *goldrefs, int numgoldrefs)
{
   int           i, j, k, l, u, v, ui, vi, x, idx;
   char          pattern[10000], *s;
   transcript_t *t1, *t2;

   for (i=start; i<end; i++) {
     t1 = transcripts+i;
     if (t1->num_exons<2) continue;

     for (j=i+1; j<=end; j++) {
       t2 = transcripts+j;
       if (t2->num_exons<2) continue;

       if (t1->num_exons+t2->num_exons<6) continue;

       k = l = 0;
       while ((k<t1->num_exons) && (l<t2->num_exons)) {
         if (t1->to[k]<t2->from[l]) k++;
         else if (t2->to[l]<t1->from[k]) l++;
         else { /* overlap; check for alternative splicing */

           if (l && k<t1->num_exons-1 && l<t2->num_exons-1 && (t1->from[k+1]<t2->to[l])) {
             for (u=1; (k+u<t1->num_exons) && (t1->from[k+u]<t2->to[l]); u++);
             u--;   /* u is the number of consecutive retained introns */

             /* PRINT */
             int from, to, type;

             /* t1 'off'; t2 'on' */
             from = t2->from[l]; to = t2->to[l];

             if ((t1->from[k]==t2->from[l]) && (t1->to[k+u]==t2->to[l]))
               type = (u>1) ? MIR_ON_TYPE : IR_ON_TYPE;
             else
               type = (u>1) ? XMIR_ON_TYPE : XIR_ON_TYPE;

             /* print 'on' form: t2 */
             s = pattern;
             sprintf(s, "%d-%d", t2->from[l], t2->to[l]);  /* match full exon(s) */
             idx = event_lookup_and_add(t2->gaid,t1->to[k],t1->from[k+u],t2->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                    t1->to[k], t1->from[k+u], Q[idx].pattern, t2->ori,
                    t2->transcript_id, t2->num_exons,
                    ((t2->ori=='+') ? (l+1) : (t2->num_exons-(l+1))), u);
             printf("%d",t2->from[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
             printf("\t%d",t2->to[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
             printf("\t%f\t%d", t2->fpkm, t2->len);

             print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

             /* now print 'off' form: t1 */

             if (t1->ori=='+') {
               from = t1->from[k]; to = t1->to[k];
             } else {
               from = t1->from[k+u]; to = t1->to[k+u];
             }

             if ((t1->from[k]==t2->from[l]) && (t1->to[k+u]==t2->to[l]))
               type = (u>1) ? MIR_OFF_TYPE : IR_OFF_TYPE;
             else
               type = (u>1) ? XMIR_OFF_TYPE : XIR_OFF_TYPE;

             s = pattern;
             sprintf(s, "%d-%d", t1->from[k], t1->to[k]);     /* match whole exon(s) pattern */
             while (*s) s++;
             for (ui=1; ui<=u; ui++) {
               sprintf(s, ",%d-%d", t1->from[k+ui], t1->to[k+ui]);
               while (*s) s++;
             }
             idx = event_lookup_and_add(t1->gaid,t1->to[k],t1->from[k+u],t1->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                    t1->to[k], t1->from[k+u], Q[idx].pattern, t1->ori,
                    t1->transcript_id, t1->num_exons,
                    ((t1->ori=='+') ? (k+1) : (t1->num_exons-(k+1))), u);
             printf("%d",t1->from[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
             printf("\t%d",t1->to[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
             printf("\t%f\t%d", t1->fpkm, t1->len);

             print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

             k += u+1; l += 1;

           } else if (k && k<t1->num_exons-1 && l<t2->num_exons-1 && (t2->from[l+1]<t1->to[k])) {
             /* t1 contains the intron */

             for (v=1; (l+v<t2->num_exons) && (t2->from[l+v]<t1->to[k]); v++);
             v--;   /* number of introns retained */


             /* PRINT */
             int from, to, type;

             /* t1 'on'; t2 'off' */
             from = t1->from[k]; to = t1->to[k];

             if ((t2->from[l]==t1->from[k]) && (t2->to[l+v]==t1->to[k]))
               type = (v>1) ? MIR_ON_TYPE : IR_ON_TYPE;
             else
               type = (v>1) ? XMIR_ON_TYPE : XIR_ON_TYPE;

             s = pattern;
             sprintf(s, "%d-%d", t1->from[k], t1->to[k]);
             idx = event_lookup_and_add(t1->gaid,t2->to[l],t2->from[l+v],t1->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                    t2->to[l], t2->from[l+v], Q[idx].pattern, t1->ori,
                    t1->transcript_id, t1->num_exons,
                    ((t1->ori=='+') ? (k+1) : (t1->num_exons-(k+1))), v);
             printf("%d",t1->from[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
             printf("\t%d",t1->to[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
             printf("\t%f\t%d", t1->fpkm, t1->len);

             print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

             /* now print 'off' form: t2 */
             if (t2->ori=='+') {
               from = t2->from[l]; to = t2->to[l];
             } else {
               from = t2->from[l+v]; to = t2->to[l+v];
             }

             if ((t2->from[l]==t1->from[k]) && (t2->to[l+v]==t1->to[k]))
               type = (v>1) ? MIR_OFF_TYPE : IR_OFF_TYPE;
             else
               type = (v>1) ? XMIR_OFF_TYPE : XIR_OFF_TYPE;

             s = pattern;
             sprintf(s, "%d-%d", t2->from[l], t2->to[l]);     /* match whole exon(s) pattern */
             while (*s) s++;
             for (vi=1; vi<=v; vi++) {
               sprintf(s, ",%d-%d", t2->from[l+vi], t2->to[l+vi]);
               while (*s) s++;
             }
             idx = event_lookup_and_add(t2->gaid,t2->to[l],t2->from[l+v],t2->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                    t2->to[l], t2->from[l+v], Q[idx].pattern, t2->ori,
                    t2->transcript_id, t2->num_exons,
                    ((t2->ori=='+') ? (l+1) : (t2->num_exons-(l+1))), v);
             printf("%d",t2->from[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
             printf("\t%d",t2->to[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
             printf("\t%f\t%d", t2->fpkm, t2->len);

             print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

             l += v+1; k += 1;

           } else {
             /* overlap, no IR */
             k++; l++;
           }
         } /* end overlap */
       } /* end criss-cross */
     } /* for j */
   } /* for i */
}

static void extract_as_ae(transcript_t *transcripts, int start, int end, tmap_t *tmaps, int numtmaps, transcript_t *goldrefs, int numgoldrefs)
{
   int           i, j, k, l, x, idx;
   char          pattern[10000], *s;
   transcript_t *t1, *t2;

   for (i=start; i<end; i++) {
     t1 = transcripts+i;
     if (t1->num_exons<3) continue;

     for (j=i+1; j<=end; j++) {
       t2 = transcripts+j;
       if (t2->num_exons<3) continue;

       k = l = 1;
       while ((k<t1->num_exons-1) && (l<t2->num_exons-1)) {
         if (t1->to[k]<t2->from[l]) k++;
         else if (t2->to[l]<t1->from[k]) l++;
         else { /* overlap; check for alternative splicing */

           if (abs(t1->from[k]-t2->from[l])<W && abs(t1->to[k]-t2->to[l])<W) {
                  ;
           } else if (overlap(t1->from[k-1],t1->to[k-1],t2->to[l-1],t2->to[l-1]) &&
                      overlap(t1->from[k+1],t1->to[k+1],t2->from[l+1],t2->to[l+1]) &&
                      ((t1->from[k]!=t2->from[l]) && (t1->to[k]!=t2->to[l])) &&
                      (max(t1->to[k-1],t2->to[l-1]) < min(t1->from[k],t2->from[l])) &&
                      (max(t1->to[k],t2->to[l])) < min(t1->from[k+1],t2->from[l+1])) {

             /* PRINT */
             /* first form */
             int from, to, type;

             from = t1->from[k]; to = t1->to[k];
             type = ((t1->to[k-1]==t2->to[l-1]) && (t1->from[k+1]==t2->from[l+1])) ?
                    AE_TYPE : XAE_TYPE;

             s = pattern;
             sprintf(s, "%d-%d", t1->from[k], t1->to[k]);  /* match full exon(s) */
             idx = event_lookup_and_add(t1->gaid,t1->from[k],t1->to[k],t1->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                    t1->from[k], t1->to[k], Q[idx].pattern, t1->ori,
                    t1->transcript_id, t1->num_exons,
                    ((t1->ori=='+') ? (k+1) : (t1->num_exons-(k+1))), 1);
             printf("%d",t1->from[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
             printf("\t%d",t1->to[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
             printf("\t%f\t%d", t1->fpkm, t1->len);

             print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

             /* now print form 2 */
             from = t2->from[l]; to = t2->to[l]; /* type stays the same */

             s = pattern;
             sprintf(s, "%d-%d", t2->from[l], t2->to[l]);  /* match full exon(s) */
             idx = event_lookup_and_add(t2->gaid,t2->from[l],t2->to[l],t2->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                    t2->from[l], t2->to[l], Q[idx].pattern, t2->ori,
                    t2->transcript_id, t2->num_exons,
                    ((t2->ori=='+') ? (l+1) : (t2->num_exons-(l+1))), 1);
             printf("%d",t2->from[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
             printf("\t%d",t2->to[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
             printf("\t%f\t%d", t2->fpkm, t2->len);

             print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

           } else if (overlap(t1->from[k-1],t1->to[k-1],t2->to[l-1],t2->to[l-1]) &&
                      (t1->to[k]==t2->to[l]) &&
                      (max(t1->to[k-1],t2->to[l-1]) < min(t1->from[k],t2->from[l])) &&
                      (max(t1->to[k],t2->to[l]) < min(t1->from[k+1],t2->from[l+1]))) {

             /* PRINT */
             /* first form */
             int from, to, type;

             from = t1->from[k]; to = t1->to[k];
             type = (t1->to[k-1]==t2->to[l-1]) ? AE_TYPE : XAE_TYPE;

             s = pattern;
             sprintf(s, "%d-%d", t1->from[k], t1->to[k]);  /* match full exon(s) */
             idx = event_lookup_and_add(t1->gaid,t1->from[k],t1->to[k],t1->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                    t1->from[k], t1->to[k], Q[idx].pattern, t1->ori,
                    t1->transcript_id, t1->num_exons,
                    ((t1->ori=='+') ? (k+1) : (t1->num_exons-(k+1))), 1);
             printf("%d",t1->from[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
             printf("\t%d",t1->to[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
             printf("\t%f\t%d", t1->fpkm, t1->len);

             print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

             /* now print form 2 */
             from = t2->from[l]; to = t2->to[l]; /* type stays the same */

             s = pattern;
             sprintf(s, "%d-%d", t2->from[l], t2->to[l]);  /* match full exon(s) */
             idx = event_lookup_and_add(t2->gaid,t2->from[l],t2->to[l],t2->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                    t2->from[l], t2->to[l], Q[idx].pattern, t2->ori,
                    t2->transcript_id, t2->num_exons,
                    ((t2->ori=='+') ? (l+1) : (t2->num_exons-(l+1))), 1);
             printf("%d",t2->from[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
             printf("\t%d",t2->to[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
             printf("\t%f\t%d", t2->fpkm, t2->len);

             print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

           } else if (overlap(t1->from[k+1],t1->to[k+1],t2->from[l+1],t2->to[l+1]) &&
                      (t1->from[k]==t2->from[l]) &&
                      (max(t1->to[k-1],t2->to[l-1]) < min(t1->from[k],t2->from[l])) &&
                      (max(t1->to[k],t2->to[l]) < min(t1->to[k+1],t2->to[l+1]))) {

             /* PRINT */
             /* first form */
             int from, to, type;

             from = t1->from[k]; to = t1->to[k];
             type = (t1->from[k+1]==t2->from[l+1]) ? AE_TYPE : XAE_TYPE; 

             s = pattern;
             sprintf(s, "%d-%d", t1->from[k], t1->to[k]);  /* match full exon(s) */
             idx = event_lookup_and_add(t1->gaid,t1->from[k],t1->to[k],t1->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t1->gene_id, Chroms[t1->gaid],
                    t1->from[k], t1->to[k], Q[idx].pattern, t1->ori,
                    t1->transcript_id, t1->num_exons,
                    ((t1->ori=='+') ? (k+1) : (t1->num_exons-(k+1))), 1);
             printf("%d",t1->from[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->from[x]);
             printf("\t%d",t1->to[0]);
             for (x=1; x<t1->num_exons; x++) printf(",%d", t1->to[x]);
             printf("\t%f\t%d", t1->fpkm, t1->len);

             print_goldref_context(t1->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);

             /* now print form 2 */
             from = t2->from[l]; to = t2->to[l]; /* type stays the same */

             s = pattern;
             sprintf(s, "%d-%d", t2->from[l], t2->to[l]);  /* match full exon(s) */
             idx = event_lookup_and_add(t2->gaid,t2->from[l],t2->to[l],t2->ori,type,pattern);

             printf("%d\t%s\t%s\t%s\t%d\t%d\t%s\t%c\t%s\t%d\t%d\t%d\t",
                    Q[idx].event_id, Codes[Q[idx].type], t2->gene_id, Chroms[t2->gaid],
                    t2->from[l], t2->to[l], Q[idx].pattern, t2->ori,
                    t2->transcript_id, t2->num_exons,
                    ((t2->ori=='+') ? (l+1) : (t2->num_exons-(l+1))), 1);
             printf("%d",t2->from[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->from[x]);
             printf("\t%d",t2->to[0]);
             for (x=1; x<t2->num_exons; x++) printf(",%d", t2->to[x]);
             printf("\t%f\t%d", t2->fpkm, t2->len);

             print_goldref_context(t2->transcript_id,from,to,tmaps,numtmaps,goldrefs,numgoldrefs);
           }
           ++k; ++l;
         } /* overlap */
       } /* end criss-cross */
     }  /* end j */
   }  /* end i */
}

/****************      Utilities      ****************/
static int lookup_chrom(char *subj, char *chroms[], int numchroms)
{
   int i;

   for (i=0; (i<numchroms) && strcmp(chroms[i],subj); i++);
   return (i==numChroms) ? -1 : i;
}

static int lookup_tmap(char *transcript_id, tmap_t *tmaps, int numtmaps)
{
   return lookup_tmap2(transcript_id, tmaps, 0, numtmaps-1);
}

static int lookup_tmap2(char *transcript_id, tmap_t *tmaps, int l, int u)
{
   int m = (int)((l+u)/2);
   int st;

   if (l>u) return -1;

   if (l==u)
     return (!strcmp(transcript_id, tmaps[l].txpt_id) ? l : -1);


   st = strcmp(transcript_id, tmaps[m].txpt_id);

   if (!st) return m;
   else if (st<0) return lookup_tmap2(transcript_id,tmaps,l,m-1);
   else return lookup_tmap2(transcript_id,tmaps,m+1,u);
}

static int lookup_goldref(char *goldref_id, transcript_t *goldrefs, int numgoldrefs)
{
   /* gaid limited list, sorted by gene_id */
   return lookup_goldref2(goldref_id, goldrefs, 0, numgoldrefs-1); 
}

static int lookup_goldref2(char *goldref_id, transcript_t *goldrefs, int l, int u)
{
   int m = (int)((l+u)/2);
   int st;

   if (l>u) return -1;

   if (l==u)
     return (!strcmp(goldref_id, goldrefs[l].transcript_id) ? l : -1);
 
   st = strcmp(goldref_id, goldrefs[m].transcript_id);
   if (!st) return m;
   if (st<0) return lookup_goldref2(goldref_id,goldrefs,l,m-1);
   else return lookup_goldref2(goldref_id,goldrefs,m+1,u);
}

static void init_event(void)
{
   nQ = 0;
}

static int  event_lookup_and_add(int gaid, int from, int to, char ori, int type, char *pattern)
{
  int i;

  for (i=0; i<nQ; i++)
    if ((Q[i].gaid == gaid) && (Q[i].ori == ori) && (Q[i].type == type) &&
        !strcmp(Q[i].pattern,pattern))
        return i;

  /* add to the list */
  if (nQ>=maxQ) {
    maxQ = max(DEFAULT_ALLOC_SMALL,2*maxQ);
    Q = realloc(Q, maxQ*sizeof(event_t));
  }

  Q[nQ].from = from;
  Q[nQ].to = to;
  Q[nQ].gaid = gaid;
  Q[nQ].type = type;
  Q[nQ].ori = ori;
  strcpy(Q[nQ].pattern,pattern);
  Q[nQ].event_id = ++eventCounter;

  nQ++;

  return (nQ-1);
}

static int alternative_tss(transcript_t *t1, transcript_t *t2, int type)
{
  /* have one internal exon edge in common and different starts */

  int i, j, found, pos;
  transcript_t *pst1, *pst2;

  if (((type != TSS_TYPE) && (type != TTS_TYPE)) || strcmp(t1->gene_id,t2->gene_id) ||
       (t1->gaid!=t2->gaid) || (t1->ori!=t2->ori) ||
       (t1->num_exons<2) || (t2->num_exons<2))
     return 0;

  /* check for shared edge */ 
  found = 0;
  i = j = 1;
  while (!found && (i<t1->num_exons) && (j<t2->num_exons)) {
    if (t1->from[i]==t2->from[j])
      found = 1;
    else if (t1->from[i]<t2->from[j])
      i++;
    else j++;
  }
  i = j = 0; 
  while (!found && (i<t1->num_exons-1) && (j<t2->num_exons-1)) {
    if (t1->to[i]==t2->to[j])
      found = 1;
    else if (t1->to[i]<t2->to[j])
      i++;
    else j++;
  }

  if (!found) return 0;

  /* shared internal junction; proceed with checking for TSS/TTS patterns */ 
  if (((type == TSS_TYPE) && (t1->ori == '+')) || ((type == TTS_TYPE) && (t1->ori == '-'))) {
    if (t1->to[0]>=t2->to[0]) {
      pst1 = t1; pst2 = t2;
    } else {
      pst1 = t2; pst2 = t1;
    }
    pos = pst1->to[0];

    j = 0; while ((j<pst2->num_exons) && (pst2->to[j]<=pos)) j++; j--;
    assert(j>=0);
      
    if (pst2->to[j]<pos) { assert (j<pst2->num_exons-1); return 1; }

    /* distinct splicing patterns and long enough overhangs */
#if 1
    assert(pst2->to[j]==pos);
    if (j && (pst2->from[j]-pst1->from[0]>V)) return 1;
#endif

  } else {
    /* t1->ori == '-' */
    if (t1->from[t1->num_exons-1]<=t2->from[t2->num_exons-1]) {
      pst1 = t1; pst2 = t2;
    } else {
      pst1 = t2; pst2 = t1;
    }
    pos = pst1->from[pst1->num_exons-1];

    j = pst2->num_exons-1; while ((j>=0) && (pst2->from[j]>=pos)) j--; j++;
    assert(j<=pst2->num_exons-1);

    if (pst2->from[j]>pos) { assert (j>0); return 1; }

    /* distinct splicing patterns and long enough overhangs */
#if 1
    assert(pst2->from[j]==pos);
    if ((j<pst2->num_exons-1) && (pst1->to[0]-pst2->to[j]>V)) return 1;
#endif
  }

  return 0;
}


/* both xi and xo are counted from 1, and starting at the 5' end */ 
static void exon_match(int from, int to, transcript_t *g, char *sign, int *xo)
{
  int i;

  if (g->ori == '+') {

    if (to<g->from[0]) {
      *sign = '<'; *xo = 0+1;
    } else if (from>g->to[g->num_exons-1]) {
      *sign = '>'; *xo = g->num_exons;
    } else {
      for (i=0; i<g->num_exons; i++) {
        if (overlap(from, to, g->from[i], g->to[i])) {
          *sign = (from==g->from[i] && to==g->to[i]) ? '=' : '~';
          *xo = i+1; 

          break;
        } else if (g->to[i]<from && to<g->from[i+1]) {
          *sign = '>';
          *xo = i+1;
 
          break;
        }
      }
    }
  } else {
    /* g->ori == '-' */

    if (to<g->from[0]) {
      *sign = '>'; *xo = g->num_exons;
    } else if (from>g->to[g->num_exons-1]) {
      *sign = '<'; *xo = 0+1;
    } else {
      for (i=g->num_exons-1; i>=0; i--) {
        if (overlap(from, to, g->from[i], g->to[i])) {
          *sign = (from==g->from[i] && to==g->to[i]) ? '=' : '~';
          *xo = g->num_exons-i;

          break;
        } else if (g->to[i]<from && to<g->from[i+1]) {
          *sign = '>';
          *xo = g->num_exons-i;
        
          break;
        }
      }
    }
  }

  return;
}

static int overlap(int a, int b, int x, int y)
{
  return ((a<=x && x<=b) || (x<=a && a<=y));
}
