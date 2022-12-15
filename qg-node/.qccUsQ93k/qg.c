#line 0 "qg-cpp.c"
#line 0 "<built-in>"
#line 0 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 0 "<command-line>"
#line 1 "qg-cpp.c"
#if _XOPEN_SOURCE < 700
  #undef _XOPEN_SOURCE
  #define _XOPEN_SOURCE 700
#endif
#if 1
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif


#line 1 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/stdlib.h"
#include <stdlib.h>
#line 2 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 4 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 5 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/stdarg.h"
#include <stdarg.h>
#line 6 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/string.h"
#include <string.h>
#line 7 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/float.h"
#include <float.h>
#line 8 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/limits.h"
#include <limits.h>
#line 9 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/math.h"
#include <math.h>
#line 10 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/time.h"
#include <time.h>
#line 11 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/sys/time.h"
#include <sys/time.h>
#line 12 "/home/lennard/basilisk/src/common.h"
#line 1 "/home/lennard/basilisk/src/ast/std/sys/resource.h"
#include <sys/resource.h>
#line 13 "/home/lennard/basilisk/src/common.h"

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif _MPI

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif
#line 46 "/home/lennard/basilisk/src/common.h"
#undef HUGE
#define HUGE ((double)1e30)

#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : HUGE)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))

#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if _MPI
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 80

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#line 104 "/home/lennard/basilisk/src/common.h"
#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if 1
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  if (!(d != NULL)) qassert ("/home/lennard/basilisk/src/common.h", 179, "d != NULL");
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", ferr);
    if (d->size == 0)
      fputs (", possible double free()", ferr);
    else
      fputs (", not traced?", ferr);
    fputs (", aborting...\n", ferr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (ferr,
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (ferr, "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (ferr,
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (ferr, "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", ferr);
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  if (!(t->stack.len > 0)) qassert ("/home/lennard/basilisk/src/common.h", 451, "t->stack.len > 0");
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define tracing(func, file, line) trace_push (&trace_func, func)
# define end_tracing(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if _MPI
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  if (!(Trace.stack.len >= 2*sizeof(double))) qassert ("/home/lennard/basilisk/src/common.h", 555, "Trace.stack.len >= 2*sizeof(double)");
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if _MPI
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if _MPI
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define tracing(...)
# define end_tracing(...)
#endif



#if _OPENMP

#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#elif _MPI

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/home/lennard/basilisk/src/common.h", 665, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 667

#define prof_stop()\
  if (!(in_prof)) qassert ("/home/lennard/basilisk/src/common.h", 669, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 672


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)
#else
     
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{tracing("mpi_all_reduce0","/home/lennard/basilisk/src/common.h",679);
  { int _ret= MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);end_tracing("mpi_all_reduce0","/home/lennard/basilisk/src/common.h",682);return _ret;}
end_tracing("mpi_all_reduce0","/home/lennard/basilisk/src/common.h",683);}
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 691

#define mpi_all_reduce_array(v,type,op,elem) {\
  prof_start ("mpi_all_reduce");\
  type global[elem], tmp[elem];\
  for (int i = 0; i < elem; i++)\
    tmp[i] = (v)[i];\
  MPI_Datatype datatype;\
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;\
  else if (!strcmp(#type, "int")) datatype = MPI_INT;\
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;\
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;\
  else {\
    fprintf (stderr, "unknown reduction type '%s'\n", #type);\
    fflush (stderr);\
    abort();\
  }\
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);\
  for (int i = 0; i < elem; i++)\
    (v)[i] = global[i];\
  prof_stop();\
}\

#line 712


#endif

#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#endif

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#line 823 "/home/lennard/basilisk/src/common.h"
#if (1 || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/lennard/basilisk/src/common.h", 833, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;


int N = 64;




typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;




} vector;

typedef struct {
  scalar * x;

  scalar * y;




} vectorl;

typedef struct {
  vector x;

  vector y;




} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 917 "/home/lennard/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  
    norm += sq(n->x);
    
#line 921
norm += sq(n->y);
  norm = sqrt(norm);
  
    n->x /= norm;
    
#line 924
n->y /= norm;
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }






  enum { right, left, top, bottom };



int nboundary = 2*2;



#define dirichlet(expr) (2.*(expr) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define dirichlet_face(expr) (expr)
#define dirichlet_face_homogeneous() (0.)
#define neumann(expr) (Delta*(expr) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#line 1 "/home/lennard/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/lennard/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 964 "/home/lennard/basilisk/src/common.h"



typedef struct {
  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

  
#line 19 "/home/lennard/basilisk/src/grid/stencils.h"
bool input, output;
  int width;
  int dirty;
  
#line 18 "/home/lennard/basilisk/src/grid/multigrid-common.h"
void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);

#line 987 "/home/lennard/basilisk/src/common.h"
} _Attributes;

static _Attributes * _attribute = NULL;






int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ ns++;}}
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s1=*_i;(&s1)->i>=0;s1=*++_i){
      if (s1.i == s.i)
 return true;}}
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      list = list_append (list, s);}}
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  {scalar*_i=(scalar*)( l2);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    l3 = list_append (l3, s);}}
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);}}
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ nv++;}}
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  {vector*_i=(vector*)( list);if(_i)for(vector w=*_i;(&w)->x.i>=0;w=*++_i){ {
    bool id = true;
    
      if (w.x.i != v.x.i)
 id = false;
      
#line 1088
if (w.y.i != v.y.i)
 id = false;
    if (id)
      return list;
  }}}
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    {vector*_i=(vector*)( l);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      list = vectors_append (list, v);}}
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
     {
      if (!(s->i >= 0)) qassert ("/home/lennard/basilisk/src/common.h", 1111, "s->i >= 0");
      v.x = *s++;
    } 
#line 1110
{
      if (!(s->i >= 0)) qassert ("/home/lennard/basilisk/src/common.h", 1111, "s->i >= 0");
      v.y = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  {tensor*_i=(tensor*)( list);if(_i)for(tensor t=*_i;(&t)->x.x.i>=0;t=*++_i){ nt++;}}
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
     {
      if (!(v->x.i >= 0)) qassert ("/home/lennard/basilisk/src/common.h", 1142, "v->x.i >= 0");
      t.x = *v++;
    } 
#line 1141
{
      if (!(v->y.i >= 0)) qassert ("/home/lennard/basilisk/src/common.h", 1142, "v->x.i >= 0");
      t.y = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
static void _init_solver (void);

void init_solver()
{
  Events = pmalloc (sizeof (Event),__func__,__FILE__,__LINE__);
  Events[0].last = 1;
  _attribute = pcalloc (datasize/sizeof(double), sizeof (_Attributes),__func__,__FILE__,__LINE__);
  int n = datasize/sizeof(double);
  all = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;
#if _CADNA
  cadna_init (-1);
#endif
#if _MPI
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}



#if _MPI
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if _MPI
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



const vector zerof = {{_NVARMAX+0},{_NVARMAX+1}};
const vector unityf = {{_NVARMAX+2},{_NVARMAX+3}};
const scalar unity = {_NVARMAX+4};
const scalar zeroc = {_NVARMAX+5};



        vector fm = {{_NVARMAX+2},{_NVARMAX+3}};
        scalar cm = {_NVARMAX+4};
#line 1269 "/home/lennard/basilisk/src/common.h"
static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 do { double __tmp = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = __tmp; } while(0);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 do { double __tmp = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = __tmp; } while(0);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}



typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

struct _display {
  const char * commands;
  bool overwrite;
};

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (struct _display p)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (p.overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(p.commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, p.commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(p.commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, p.commands);
  }
}







#line 1 "/home/lennard/basilisk/src/grid/stencils.h"
#line 17 "/home/lennard/basilisk/src/grid/stencils.h"










typedef struct {
  const char * fname;
  int line;
  int first;
  int face;
  bool vertex;
} ForeachData;


#define foreach_stencil() {\
  static ForeachData _loop = {\
    __FILE__, __LINE__,\
    1, 0, 0\
  };\
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock;\
  s.i >= 0; i++, s = *i) {\
    _attribute[s.i].input = _attribute[s.i].output = false;\
    _attribute[s.i].width = 0;\
  }\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0}; NOT_UNUSED (point);\

#line 48


#define end_foreach_stencil()\
  end_stencil (&_loop);\
  _loop.first = 0;\
}\

#line 54


#define foreach_vertex_stencil() foreach_stencil() _loop.vertex = true;
#define end_foreach_vertex_stencil() end_foreach_stencil()

#define foreach_face_stencil() foreach_stencil()
#define end_foreach_face_stencil() end_foreach_stencil()

#define foreach_visible_stencil(...) foreach_stencil()
#define end_foreach_visible_stencil(...) end_foreach_stencil()

#define _stencil_is_face_x() { _loop.face |= (1 << 0);
#define end__stencil_is_face_x() }
#define _stencil_is_face_y() { _loop.face |= (1 << 1);
#define end__stencil_is_face_y() }
#define _stencil_is_face_z() { _loop.face |= (1 << 2);
#define end__stencil_is_face_z() }

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow);
void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line);

#define _stencil_val(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, false)\

#line 79

#define _stencil_val_o(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, true)\

#line 82

#define _stencil_val_a(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, false, __FILE__, __LINE__)\

#line 85

#define _stencil_val_r(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, true, __FILE__, __LINE__)\

#line 88


#define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_fine_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_fine_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define _stencil_coarse(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_coarse_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_coarse_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define r_assign(x)
#define _assign(x)

#define _stencil_neighbor(i,j,k)
#define _stencil_child(i,j,k)
#define _stencil_aparent(i,j,k)
#define _stencil_aparent_a(i,j,k)
#define _stencil_aparent_r(i,j,k)

#define _stencil_neighborp(i,j,k) neighborp(i,j,k)

int _stencil_nop;
#define _stencil_val_higher_dimension (_stencil_nop = 1)
#define _stencil__val_constant(a,_i,_j,_k) (_stencil_nop = 1)

typedef void _stencil_undefined;

#define o_stencil -2







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      return true;}}
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (s.i == b.i)
      return true;}}
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void end_stencil (ForeachData * loop)
{
  scalar * listc = NULL, * dirty = NULL;
  vectorl listf = {NULL};
  bool flux = false;




  {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    bool write = _attribute[s.i].output, read = _attribute[s.i].input;




    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (_attribute[s.i].width > 0)
     listc = list_append (listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     
       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.x = list_append (listf.x, s);
    flux = true;
  }
       }
       
#line 194
if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.y = list_append (listf.y, s);
    flux = true;
  }
       }
   }
 }





 else if (_attribute[s.i].width > 0)
   listc = list_append (listc, s);
      }





      if (write) {
 if (2 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 225
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
      {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     } 
#line 243
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     }
     if (!_attribute[s.i].face && loop->first)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     struct { int x, y, z; } input, output;
     vector v = _attribute[s.i].v;

     
       input.x = _attribute[v.x.i].input, output.x = _attribute[v.x.i].output;
       
#line 273
input.y = _attribute[v.y.i].input, output.y = _attribute[v.y.i].output;

     init_face_vector (v, name);


     
       _attribute[v.x.i].input = input.x, _attribute[v.x.i].output = output.x;
       
#line 279
_attribute[v.y.i].input = input.y, _attribute[v.y.i].output = output.y;





     pfree (name,__func__,__FILE__,__LINE__);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 291
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     init_vertex_scalar (s, name);
     
       _attribute[s.i].v.x.i = -1;
       
#line 298
_attribute[s.i].v.y.i = -1;




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }





 dirty = list_append (dirty, s);
 {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
   if (scalar_depends_from (d, s))
     dirty = list_append (dirty, d);}}
      }
    }
  }}}




  if (flux) {
#line 335 "/home/lennard/basilisk/src/grid/stencils.h"
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,__LINE__);
      
#line 337
pfree (listf.y,__func__,__FILE__,__LINE__);
  }




  if (listc) {






    boundary_internal (listc, loop->fname, loop->line);
    pfree (listc,__func__,__FILE__,__LINE__);
  }





  if (dirty) {






    {scalar*_i=(scalar*)( dirty);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = true;}}
    pfree (dirty,__func__,__FILE__,__LINE__);
  }
}
#line 1445 "/home/lennard/basilisk/src/common.h"
#line 13 "qg-cpp.c"
#line 1 "qg.c"
#line 20 "qg.c"
#line 1 "grid/multigrid.h"
#line 1 "/home/lennard/basilisk/src/grid/multigrid.h"







# define BGHOSTS 1







typedef struct {
  Grid g;
  char ** d;
} Multigrid;

struct _Point {
  int i;

  int j;




  int level, n;
#ifdef foreach_block
  int l;
  #define _BLOCK_INDEX , point.l
#else
  #define _BLOCK_INDEX
#endif
};
static Point last_point;
#line 49 "/home/lennard/basilisk/src/grid/multigrid.h"
static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*2;
  return sq(n);
}
#line 62 "/home/lennard/basilisk/src/grid/multigrid.h"
#define data(k,l,m)\
  ((double *)&((Multigrid *)grid)->d[point.level][((point.i + k)*((1 << point.level) +\
       2*2) +\
      (point.j + l))*datasize]) 
#line 64

#line 89 "/home/lennard/basilisk/src/grid/multigrid.h"
#define allocated(k,l,m) (point.i+k >= 0 && point.i+k < (1 << point.level) + 2*2 &&\
         point.j+l >= 0 && point.j+l < (1 << point.level) + 2*2)\

#line 91


#define allocated_child(k,l,m) (level < depth() &&\
         point.i > 0 && point.i <= (1 << point.level) + 2 &&\
         point.j > 0 && point.j <= (1 << point.level) + 2)\

#line 96

#line 117 "/home/lennard/basilisk/src/grid/multigrid.h"
#define depth() (grid->depth)
#line 136 "/home/lennard/basilisk/src/grid/multigrid.h"
#define fine(a,k,l,m)\
  ((double *)\
   &((Multigrid *)grid)->d[point.level+1][((2*point.i-2 +k)*2*((1 << point.level) +\
        2) +\
     (2*point.j-2 +l))*datasize])[_index(a,m)]\

#line 141

#define coarse(a,k,l,m)\
  ((double *)\
   &((Multigrid *)grid)->d[point.level-1][(((point.i+2)/2+k)*((1 << point.level)/2 +\
        2*2) +\
     (point.j+2)/2+l)*datasize])[_index(a,m)]\

#line 147

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
  struct { int x, y; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1\
  }; NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2; parent.j = (point.j + 2)/2;\

#line 157

#line 191 "/home/lennard/basilisk/src/grid/multigrid.h"
#define foreach_level(l)\
OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
\
\
 {\
\
          POINT_VARIABLES\

#line 208

#define end_foreach_level()\
\
 }\
\
  }\
}\

#line 215


#define foreach()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
\
\
 {\
\
          POINT_VARIABLES\

#line 234

#define end_foreach()\
\
 }\
\
  }\
}\

#line 241


#define is_active(cell) (true)
#define is_leaf(cell) (level == depth())
#define is_local(cell) (true)
#define leaf 2
#define refine_cell(...) do {\
  fprintf (stderr, "grid depths do not match. Aborting.\n");\
  if (!(0)) qassert ("/home/lennard/basilisk/src/grid/multigrid.h", 249, "0");\
} while (0)\

#line 251

#define tree ((Multigrid *)grid)
#line 1 "grid/foreach_cell.h"
#line 1 "/home/lennard/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/lennard/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
  Point root = {2,2,0};\
\
\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = {0};\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
\
\
\
\
\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
    Point root = {2,2,0};\
\
\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = {0};\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 254 "/home/lennard/basilisk/src/grid/multigrid.h"

#define foreach_face_generic()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k <= point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j <= point.n + 2; point.j++)\
\
\
\
        {\
\
   POINT_VARIABLES\

#line 272

#define end_foreach_face_generic()\
\
 }\
\
  }\
}\

#line 279


#define foreach_vertex()\
foreach_face_generic() {\
  x -= Delta/2.;\
\
  y -= Delta/2.;\
\
\
\
\

#line 290

#define end_foreach_vertex() } end_foreach_face_generic()

#define is_coarse() (point.level < depth())
#line 321 "/home/lennard/basilisk/src/grid/multigrid.h"
#define is_face_x() { int ig = -1; VARIABLES; if (point.j < point.n + 2) {
#define end_is_face_x() }}
#define is_face_y() { int jg = -1; VARIABLES; if (point.i < point.n + 2) {
#define end_is_face_y() }}

#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2;\
  point.level++;\
  point.n *= 2;\
  for (int _k = 0; _k < 2; _k++)\
    for (int _l = 0; _l < 2; _l++) {\
      point.i = _i + _k; point.j = _j + _l;\
      POINT_VARIABLES;\

#line 334

#define end_foreach_child()\
  }\
  point.i = (_i + 2)/2; point.j = (_j + 2)/2;\
  point.level--;\
  point.n /= 2;\
}\

#line 341

#define foreach_child_break() _k = _l = 2
#line 387 "/home/lennard/basilisk/src/grid/multigrid.h"
#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

#line 1 "grid/neighbors.h"
#line 1 "/home/lennard/basilisk/src/grid/neighbors.h"
#line 17 "/home/lennard/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j;\
  for (int _k = - _nn; _k <= _nn; _k++) {\
    point.i = _i + _k;\
    for (int _l = - _nn; _l <= _nn; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 25

#define end_foreach_neighbor()\
    }\
  }\
  point.i = _i; point.j = _j;\
}\

#line 31

#define foreach_neighbor_break() _k = _l = _nn + 1
#line 393 "/home/lennard/basilisk/src/grid/multigrid.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Point p;
  p.level = depth(); p.n = 1 << p.level;
  for (; p.level >= 0; p.n /= 2, p.level--)
    for (int i = 0; i < sq(p.n + 2*2); i++)
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 if (!is_constant(s))
   for (int b = 0; b < _attribute[s.i].block; b++)
     ((double *)(&((Multigrid *)grid)->d[p.level][i*datasize]))[s.i + b] = val;
      }}}
}
#line 433 "/home/lennard/basilisk/src/grid/multigrid.h"
#define foreach_boundary_dir(l,d)\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l < 0 ? depth() : l;\
  point.n = 1 << point.level;\
  int * _i = &point.j;\
  if (d == left) {\
    point.i = 2;\
    ig = -1;\
  }\
  else if (d == right) {\
    point.i = point.n + 2 - 1;\
    ig = 1;\
  }\
  else if (d == bottom) {\
    point.j = 2;\
    _i = &point.i;\
    jg = -1;\
  }\
  else if (d == top) {\
    point.j = point.n + 2 - 1;\
    _i = &point.i;\
    jg = 1;\
  }\
  int _l;\
  OMP(omp for schedule(static))\
  for (_l = 0; _l < point.n + 2*2; _l++) {\
    *_i = _l;\
    {\
      POINT_VARIABLES\

#line 464

#define end_foreach_boundary_dir()\
    }\
  }\
}\

#line 469


#define neighbor(o,p,q)\
  ((Point){point.i+o, point.j+p, point.level, point.n _BLOCK_INDEX})\

#line 473

#define is_boundary(point) (point.i < 2 || point.i >= point.n + 2 ||\
    point.j < 2 || point.j >= point.n + 2)\

#line 476

#line 538 "/home/lennard/basilisk/src/grid/multigrid.h"
#define foreach_boundary(b)\
  if (default_scalar_bc[b] != periodic_bc)\
    foreach_boundary_dir (depth(), b)\
      if (!is_boundary(point)) {\

#line 542

#define end_foreach_boundary() } end_foreach_boundary_dir()

#define neighborp(k,l,o) neighbor(k,l,o)

static double periodic_bc (Point point, Point neighbor, scalar s, void * data);

static inline bool is_vertex_scalar (scalar s)
{
  
    if (_attribute[s.i].d.x != -1)
      return false;
    
#line 552
if (_attribute[s.i].d.y != -1)
      return false;
  return true;
}

static void box_boundary_level (const Boundary * b, scalar * scalars, int l)
{
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  for (int bghost = 1; bghost <= BGHOSTS; bghost++)
    for (int d = 0; d < 2*2; d++) {

      scalar * list = NULL, * listb = NULL;
      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   scalar sb = s;

   if (_attribute[s.i].v.x.i >= 0) {

     int j = 0;
     while ((&_attribute[s.i].v.x)[j].i != s.i) j++;
     sb = (&_attribute[s.i].v.x)[(j - d/2 + 2) % 2];
   }

   if (_attribute[sb.i].boundary[d] && _attribute[sb.i].boundary[d] != periodic_bc) {
     list = list_append (list, s);
     listb = list_append (listb, sb);
   }
 }}}

      if (list) {
 extern double (* default_scalar_bc[]) (Point, Point, scalar, void *);
 if (default_scalar_bc[d] != periodic_bc)
 {foreach_boundary_dir (l, d) {
   scalar s, sb;
   {scalar*_i0=listb;scalar*_i1= list;if(_i0)for(sb=*_i0,s=*_i1;_i0->i>= 0;sb=*++_i0,s=*++_i1){ {
     if ((_attribute[s.i].face && sb.i == _attribute[s.i].v.x.i) || is_vertex_scalar (s)) {

       if (bghost == 1)
 
    val(s,(ig + 1)/2,(jg + 1)/2,(kg + 1)/2) =
    _attribute[sb.i].boundary[d] (point, neighborp(ig,jg,kg), s, NULL);
     }
     else

      
  val(s,bghost*ig,bghost*jg,bghost*kg) =
  _attribute[sb.i].boundary[d] (neighborp((1 - bghost)*ig,
       (1 - bghost)*jg,
       (1 - bghost)*kg),
    neighborp(bghost*ig,bghost*jg,bghost*kg),
    s, NULL);
   }}}
 }end_foreach_boundary_dir();}
 pfree (list,__func__,__FILE__,__LINE__);
 pfree (listb,__func__,__FILE__,__LINE__);
      }
    }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#line 652 "/home/lennard/basilisk/src/grid/multigrid.h"
#define VT _attribute[s.i].v.y


static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].face) {
 scalar vt = VT;
 if (_attribute[vt.i].boundary[right] == periodic_bc)
   list1 = list_add (list1, s);
      }
      else if (_attribute[s.i].boundary[right] == periodic_bc)
 list1 = list_add (list1, s);
    }}}
  if (!list1)
    return;

  if (l == 0) {
    {foreach_level(0)
      {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double * v = &val(s,0,0,0);
 {foreach_neighbor()
   memcpy (&val(s,0,0,0), v, _attribute[s.i].block*sizeof(double));end_foreach_neighbor()}
      }}}end_foreach_level();}
    pfree (list1,__func__,__FILE__,__LINE__);
    return;
  }

  OMP_PARALLEL() {
    Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    point.level = l < 0 ? depth() : l; point.n = 1 << point.level;

    int j;
    OMP(omp for schedule(static))
      for (j = 0; j < point.n + 2*2; j++) {
 for (int i = 0; i < 2; i++)
   {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
     memcpy (&val(s,i,j,0), &val(s,i + point.n,j,0), _attribute[s.i].block*sizeof(double));}}
 for (int i = point.n + 2; i < point.n + 2*2; i++)
   {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
     memcpy (&val(s,i,j,0), &val(s,i - point.n,j,0), _attribute[s.i].block*sizeof(double));}}
      }
#line 709 "/home/lennard/basilisk/src/grid/multigrid.h"
  }
  pfree (list1,__func__,__FILE__,__LINE__);
}

#line 655
static void periodic_boundary_level_y (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].face) {
 scalar vt = VT;
 if (_attribute[vt.i].boundary[top] == periodic_bc)
   list1 = list_add (list1, s);
      }
      else if (_attribute[s.i].boundary[top] == periodic_bc)
 list1 = list_add (list1, s);
    }}}
  if (!list1)
    return;

  if (l == 0) {
    {foreach_level(0)
      {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double * v = &val(s,0,0,0);
 {foreach_neighbor()
   memcpy (&val(s,0,0,0), v, _attribute[s.i].block*sizeof(double));end_foreach_neighbor()}
      }}}end_foreach_level();}
    pfree (list1,__func__,__FILE__,__LINE__);
    return;
  }

  OMP_PARALLEL() {
    Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    point.level = l < 0 ? depth() : l; point.n = 1 << point.level;

    int j;
    OMP(omp for schedule(static))
      for (j = 0; j < point.n + 2*2; j++) {
 for (int i = 0; i < 2; i++)
   {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
     memcpy (&val(s,j,i,0), &val(s,j,i + point.n,0), _attribute[s.i].block*sizeof(double));}}
 for (int i = point.n + 2; i < point.n + 2*2; i++)
   {scalar*_i=(scalar*)( list1);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
     memcpy (&val(s,j,i,0), &val(s,j,i - point.n,0), _attribute[s.i].block*sizeof(double));}}
      }
#line 709 "/home/lennard/basilisk/src/grid/multigrid.h"
  }
  pfree (list1,__func__,__FILE__,__LINE__);
}

#undef VT





void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Multigrid * m = ((Multigrid *)grid);
  for (int l = 0; l <= depth(); l++)
    pfree (m->d[l],__func__,__FILE__,__LINE__);
  pfree (m->d,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
  grid = NULL;
}

int log_base2 (int n) {
  int m = n, r = 0;
  while (m > 1)
    m /= 2, r++;
  return (1 << r) < n ? r + 1 : r;
}

void init_grid (int n)
{
  free_grid();
  Multigrid * m = ((Multigrid *) pmalloc ((1)*sizeof(Multigrid),__func__,__FILE__,__LINE__));
  grid = (Grid *) m;
  grid->depth = grid->maxdepth = log_base2(n);
  N = 1 << depth();

  grid->n = grid->tn = 1 << 2*depth();

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  add_boundary (b);





   {
    Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
    b->level = periodic_boundary_level_x;
    add_boundary (b);
  } 
#line 757
{
    Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
    b->level = periodic_boundary_level_y;
    add_boundary (b);
  }


  m->d = (char **) pmalloc(sizeof(Point *)*(depth() + 1),__func__,__FILE__,__LINE__);
  for (int l = 0; l <= depth(); l++) {
    size_t len = _size(l)*datasize;
    m->d[l] = (char *) pmalloc (len,__func__,__FILE__,__LINE__);


    double * v = (double *) m->d[l];
    for (int i = 0; i < len/sizeof(double); i++)
      v[i] = undefined;
  }
  reset (all, 0.);
}

void realloc_scalar (int size)
{
  Multigrid * p = ((Multigrid *)grid);
  for (int l = 0; l <= depth(); l++) {
    size_t len = _size(l);
    p->d[l] = (char *) prealloc (p->d[l], (len*(datasize + size))*sizeof(char),__func__,__FILE__,__LINE__);
    char * data = p->d[l] + (len - 1)*datasize;
    for (int i = len - 1; i > 0; i--, data -= datasize)
      memmove (data + i*size, data, datasize);
  }
  datasize += size;
}
#line 802 "/home/lennard/basilisk/src/grid/multigrid.h"
struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  point.level = -1, point.n = 1 << depth();
#line 823 "/home/lennard/basilisk/src/grid/multigrid.h"
  point.i = (p.x - X0)/L0*point.n + 2;
  if (point.i < 2 || point.i >= point.n + 2)
    return point;

  point.j = (p.y - Y0)/L0*point.n + 2;
  if (point.j < 2 || point.j >= point.n + 2)
    return point;







  point.level = depth();
  return point;
}

#line 1 "grid/multigrid-common.h"
#line 1 "/home/lennard/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/lennard/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/lennard/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/lennard/basilisk/src/grid/events.h", 87, "Events");
  if (!(!event.last)) qassert ("/home/lennard/basilisk/src/grid/events.h", 88, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/lennard/basilisk/src/grid/events.h", 92, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/home/lennard/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/lennard/basilisk/src/grid/events.h", 245, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/lennard/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
#define diagonalize(a)
#define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  double Delta_x = Delta;\
\
  double Delta_y = Delta;\
\
\
\
\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
\
\
  double z = 0.;\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  NOT_UNUSED(Delta_x);\
\
  NOT_UNUSED(Delta_y);\
\
\
\
\
\
  ;\

#line 44


#line 1 "grid/fpe.h"
#line 1 "/home/lennard/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 47 "/home/lennard/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()



static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    _attribute[sb.i].block = block;
    init_scalar (sb, bname);
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
    init_scalar (sb, bname);
  }
  all = list_append (all, sb);
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
 init_block_scalar (sb, name, ext, n, block);
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/lennard/basilisk/src/grid/cartesian-common.h", 91, "nvar + block <= _NVARMAX");
  _attribute = (_Attributes *) prealloc (_attribute, (nvar + block)*sizeof(_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_scalar (const char * name)
{
  return new_block_scalar (name, "", 1);
}

scalar new_vertex_scalar (const char * name)
{
  scalar s = new_block_scalar (name, "", 1);
  init_vertex_scalar (s, NULL);
  return s;
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  
    v.x = new_block_scalar (name, ext.x, block);
    
#line 121
v.y = new_block_scalar (name, ext.y, block);
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 145
vb.y.i = v.y.i + i;
    init_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 148
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 151
_attribute[v.y.i].block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 161
vb.y.i = v.y.i + i;
    init_face_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 164
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 167
_attribute[v.y.i].block = block;
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  } 
#line 176
{
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  } 
#line 189
{
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;
#line 209 "/home/lennard/basilisk/src/grid/cartesian-common.h"
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  
    init_const_scalar (v.x, name, *val++);
    
#line 234
init_const_scalar (v.y, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  
    v.x.i = _NVARMAX + i++;
    
#line 241
v.y.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[a.i].boundary_homogeneous;
  if (!(_attribute[b.i].block > 0 && _attribute[a.i].block == _attribute[b.i].block)) qassert ("/home/lennard/basilisk/src/grid/cartesian-common.h", 252, "b.block > 0 && a.block == b.block");
  pfree (_attribute[a.i].depends,__func__,__FILE__,__LINE__);
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
  _attribute[a.i].depends = list_copy (_attribute[b.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) :
      new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }}}
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    {
      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
      
#line 280
if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];}}}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,__LINE__); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }}}

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    if (_attribute[f.i].block > 0) {
      scalar * s = all;
      for (; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }}}
}

void free_solver()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/home/lennard/basilisk/src/grid/cartesian-common.h", 331, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (baseblock,__func__,__FILE__,__LINE__); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    {
      list1.x = list_append (list1.x, v.x);
      
#line 383
list1.y = list_append (list1.y, v.y);}}}
  boundary_face (list1);
  
    pfree (list1.x,__func__,__FILE__,__LINE__);
    
#line 386
pfree (list1.y,__func__,__FILE__,__LINE__);
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  scalar * list1 = list;
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);}}
  return list_append (list1, s);
}

     
void boundary_internal (scalar * list, const char * fname, int line)
{tracing("boundary_internal","/home/lennard/basilisk/src/grid/cartesian-common.h",402);
  if (list == NULL)
    {end_tracing("boundary_internal","/home/lennard/basilisk/src/grid/cartesian-common.h",405);return;}
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;
     
#line 414
if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }}}
  if (flux) {
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,__LINE__);
      
#line 429
pfree (listf.y,__func__,__FILE__,__LINE__);
  }
  if (listc) {
    boundary_level (listc, -1);
    {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = false;}}
    pfree (listc,__func__,__FILE__,__LINE__);
  }
end_tracing("boundary_internal","/home/lennard/basilisk/src/grid/cartesian-common.h",437);}

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  
    {scalar*_i=(scalar*)( list.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 447
{scalar*_i=(scalar*)( list.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[s.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[s.i].boundary_homogeneous;

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;
  _attribute[s.i].block = block == 0 ? 1 : block;

  _attribute[s.i].boundary = boundary ? boundary :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*2 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
   {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  } 
#line 494
{
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  
    _attribute[s.i].d.x = -1;
    
#line 506
_attribute[s.i].d.y = -1;
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar, void *) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  } 
#line 521
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*2 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
   {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  } 
#line 541
{
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  } 
#line 553
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }






    for (int b = 0; b < nboundary; b++) {
      _attribute[t.x.x.i].boundary[b] = _attribute[t.y.x.i].boundary[b] =
 _attribute[t.x.x.i].boundary_homogeneous[b] = _attribute[t.y.y.i].boundary_homogeneous[b] =
 b < 2*2 ? default_scalar_bc[b] : symmetry;
      _attribute[t.x.y.i].boundary[b] = _attribute[t.y.y.i].boundary[b] =
 _attribute[t.x.y.i].boundary_homogeneous[b] = _attribute[t.y.x.i].boundary_homogeneous[b] =
 b < 2*2 ? default_vector_bc[b] : antisymmetry;
    }



  return t;
}

struct OutputCells {
  FILE * fp;
  coord c;
  double size;
};

void output_cells (struct OutputCells p)
{
  if (!p.fp) p.fp = fout;
  {foreach() {
    bool inside = true;
    coord o = {x,y,z};
    
      if (inside && p.size > 0. &&
   (o.x > p.c.x + p.size || o.x < p.c.x - p.size))
 inside = false;
      
#line 595
if (inside && p.size > 0. &&
   (o.y > p.c.y + p.size || o.y < p.c.y - p.size))
 inside = false;
    if (inside) {
      Delta /= 2.;



      fprintf (p.fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 623 "/home/lennard/basilisk/src/grid/cartesian-common.h"
    }
  }end_foreach();}
  fflush (p.fp);
}
#line 635 "/home/lennard/basilisk/src/grid/cartesian-common.h"
static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"




    "  plot '%s' w l lc 0, "
    "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",





    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells ((struct OutputCells){fp, (coord){x,y,z}, 4.*Delta});
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){



    fprintf (fp, "x y %s ", _attribute[v.i].name);}}



  fputc ('\n', fp);
#line 702 "/home/lennard/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
   fprintf (fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   if (allocated(k,l,0))
     fprintf (fp, "%g ", val(v,k,l,0));
   else
     fputs ("n/a ", fp);
 }}}
 fputc ('\n', fp);
      }
#line 732 "/home/lennard/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }}}
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  debug = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  scalar v = p.v;







  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);
#line 806 "/home/lennard/basilisk/src/grid/cartesian-common.h"
}

     
double interpolate (struct _interpolate p)
{tracing("interpolate","/home/lennard/basilisk/src/grid/cartesian-common.h",809);
  scalar v = p.v;
  boundary_internal ((scalar *)((scalar[]){v,{-1}}), "/home/lennard/basilisk/src/grid/cartesian-common.h", 812);
  Point point = locate ((struct _locate){p.x, p.y, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (point.level < 0)
    {end_tracing("interpolate","/home/lennard/basilisk/src/grid/cartesian-common.h",815);return HUGE;}
  { double _ret= interpolate_linear (point, p);end_tracing("interpolate","/home/lennard/basilisk/src/grid/cartesian-common.h",816);return _ret;}
end_tracing("interpolate","/home/lennard/basilisk/src/grid/cartesian-common.h",817);}

     
void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{tracing("interpolate_array","/home/lennard/basilisk/src/grid/cartesian-common.h",820);
  boundary_internal ((scalar *)list, "/home/lennard/basilisk/src/grid/cartesian-common.h", 822);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    if (point.level >= 0) {
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});}}
    }
    else
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = HUGE;}}
  }
#if _MPI
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
end_tracing("interpolate_array","/home/lennard/basilisk/src/grid/cartesian-common.h",844);}



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }}}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      
 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
 
#line 865
_attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }}}
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return HUGE;
}

static void periodic_boundary (int d)
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;}}

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }}}

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{



    if (!(dir <= bottom)) qassert ("/home/lennard/basilisk/src/grid/cartesian-common.h", 901, "dir <= bottom");




  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,i,j,k);
}

void default_stencil (Point p, scalar * list)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].input = true, _attribute[s.i].width = 2;}}
}




static void write_stencil_index (int * index)
{
  fprintf (qstderr(), "[%d", index[0]);
  for (int d = 1; d < 2; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  bool central = true;
  for (int d = 0; d < 2; d++) {
    if (!overflow && (index[d] > 2 || index[d] < - 2)) {
      fprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
    if (index[d] != 0)
      central = false;
  }
  if (central) {
    if (!_attribute[s.i].output)
      _attribute[s.i].input = true;
  }
  else {
    _attribute[s.i].input = true;
    int d = 0;
     {
      if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 963
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    }
  }
}

void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    abort();
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  for (int d = 0; d < 2; d++)
    if (index[d] != 0) {
      fprintf (qstderr(), "%s:%d: error: illegal write: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
  if (input && !_attribute[s.i].output)
    _attribute[s.i].input = true;
  _attribute[s.i].output = true;
}
#line 4 "/home/lennard/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2);
}

static inline void restriction_volume_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 35
if(!is_constant(cm)){{
  double sum = 0.;
  {foreach_child()
    sum += val(cm,0,0,0)*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(val(cm,0,0,0) + 1e-30);
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 35
{
  double sum = 0.;
  {foreach_child()
    sum += _const_cm*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(_const_cm + 1e-30);
}}

#line 40
}

static inline void face_average (Point point, vector v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
   {




      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0))/2.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0))/2.;






  } 
#line 44
{




      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,1,0,0))/2.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,1,2,0))/2.;






  }
}

static inline void restriction_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);





  }
}

static inline void no_restriction (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;}

static inline void no_data (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = HUGE;end_foreach_child()}
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar[]){s,{-1}}));
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    {foreach_coarse_level (l) {
      {foreach_child()
        val(w,0,0,0) = val(s,0,0,0);end_foreach_child()}
      _attribute[s.i].prolongation (point, s);
      {foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      }end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){w,{-1}}), l + 1);
  }

  {foreach_level(0)
    val(w,0,0,0) = val(s,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
  {foreach_level(0)
    val(s,0,0,0) = val(w,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){s,{-1}}), 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    {foreach_coarse_level (l) {
      _attribute[s.i].prolongation (point, s);
      {foreach_child()
        val(s,0,0,0) += val(w,0,0,0);end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



    return (9.*coarse(s,0,0,0) +
     3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
     coarse(s,child.x,child.y,0))/16.;
#line 140 "/home/lennard/basilisk/src/grid/multigrid-common.h"
}

static inline void refine_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = bilinear (point, s);end_foreach_child()}
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return
    quadratic (quadratic (coarse(s,0,0,0),
     coarse(s,child.x,0,0),
     coarse(s,-child.x,0,0)),
        quadratic (coarse(s,0,child.y,0),
     coarse(s,child.x,child.y,0),
     coarse(s,-child.x,child.y,0)),
        quadratic (coarse(s,0,-child.y,0),
     coarse(s,child.x,-child.y,0),
     coarse(s,-child.x,-child.y,0)));




}

static inline double biquadratic_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return (36.*val(s,0,0,0) + 18.*(val(s,-1,0,0) + val(s,0,-1,0)) - 6.*(val(s,1,0,0) + val(s,0,1,0)) +
   9.*val(s,-1,-1,0) - 3.*(val(s,1,-1,0) + val(s,-1,1,0)) + val(s,1,1,0))/64.;




}

static inline void refine_biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = biquadratic (point, s);end_foreach_child()}
}

static inline void refine_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 194
if(!is_constant(cm)){{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val(cm,0,0,0), sum = val(cm,0,0,0)*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;
    sum -= val(cm,0,0,0);
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/lennard/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 194
{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*_const_cm, sum = _const_cm*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*_const_cm/cmc;
    sum -= _const_cm;
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/lennard/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}

#line 211
}

static inline void refine_reset (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(v,0,0,0) = 0.;end_foreach_child()}
}

static inline void refine_injection (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double val = val(v,0,0,0);
  {foreach_child()
    val(v,0,0,0) = val;end_foreach_child()}
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.y.i].restriction = no_restriction;
    
#line 245
_attribute[v.x.i].restriction = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 271 "/home/lennard/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){
     fprintf (fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       coarse(v,k*child.x,l*child.y,0));}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 302 "/home/lennard/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 324 "/home/lennard/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
     fprintf (fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     if (allocated_child(k,l,0))
       fprintf (fp, "%g ", fine(v,k,l,0));
     else
       fputs ("n/a ", fp);
   }}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 362 "/home/lennard/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 381
list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }}}

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
      {foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
  
     _attribute[s.i].restriction (point, s);
 }}}
      }end_foreach_coarse_level();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




  foreach_stencil()
    {_stencil_val_a(size,0,0,0);  }end_foreach_stencil();




  {
#line 428
foreach()
    val(size,0,0,0) = 1;end_foreach();}





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
    {foreach_coarse_level(l) {
      double sum = !leaves;
      {foreach_child()
 sum += val(size,0,0,0);end_foreach_child()}
      val(size,0,0,0) = sum;
    }end_foreach_coarse_level();}
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){size,{-1}}), l); };
  }
}
#line 842 "/home/lennard/basilisk/src/grid/multigrid.h"

struct Dimensions {
  int nx, ny, nz;
};

void dimensions (struct Dimensions p)
{




}
#line 21 "qg.c"
#line 1 "qg.h"
#line 1 "./qg.h"
#line 77 "./qg.h"
#line 1 "predictor-corrector.h"
#line 1 "/home/lennard/basilisk/src/predictor-corrector.h"


#line 1 "utils.h"
#line 1 "/home/lennard/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if _MPI
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:n)){
#line 69
foreach() n++;end_foreach();mpi_all_reduce_array(&n,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
    
#line 70
s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if 1
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if _MPI
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Multigrid"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if _MPI
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach_stencil(
)
    {_stencil_val(f,0,0,0);_stencil_val(cm,0,0,0); {   
      _stencil_val(f,0,0,0);   
         
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       
    
#line 143
}       }end_foreach_stencil();
  
#line 135
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != HUGE && (sq(Delta)*val(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val(cm,0,0,0));
      avg += (sq(Delta)*val(cm,0,0,0))*v;
      rms += (sq(Delta)*val(cm,0,0,0))*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != HUGE && (sq(Delta)*_const_cm) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*_const_cm);
      avg += (sq(Delta)*_const_cm)*v;
      rms += (sq(Delta)*_const_cm)*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach_stencil(
)
    {_stencil_val(cm,0,0,0); _stencil_val(f,0,0,0); {
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
  
#line 163
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != HUGE) {
      volume += (sq(Delta)*val(cm,0,0,0));
      sum += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*_const_cm) > 0. && val(f,0,0,0) != HUGE) {
      volume += (sq(Delta)*_const_cm);
      sum += (sq(Delta)*_const_cm)*val(f,0,0,0);
      sum2 += (sq(Delta)*_const_cm)*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 187 "/home/lennard/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 213 "/home/lennard/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 237 "/home/lennard/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/home/lennard/basilisk/src/utils.h", 239, "list_len(f) == vectors_len(g)");
  foreach_stencil() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





     _stencil_val_a(v.x,0,0,0);_stencil_val(s,-1,0,0); _stencil_val(s,0,0,0); _stencil_val(s,1,0,0);   
 } 
#line 244
{





     _stencil_val_a(v.y,0,0,0);_stencil_val(s,0,-1,0); _stencil_val(s,0,0,0); _stencil_val(s,0,1,0);   
 }}
      else
 { {





     _stencil_val_a(v.x,0,0,0);_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);   
 } 
#line 253
{





     _stencil_val_a(v.y,0,0,0);_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);   
 }}
    }}}
  }end_foreach_stencil();
  {
#line 240
foreach() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





     val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
 } 
#line 244
{





     val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
 }}
      else
 { {





     val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
 } 
#line 253
{





     val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
 }}
    }}}
  }end_foreach();}
}
#line 280 "/home/lennard/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
  foreach_stencil()
    {_stencil_val_a(omega,0,0,0);_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);      
             
#line 286
}end_foreach_stencil();
  
#line 282
if(!is_constant(fm.x) && !is_constant(cm)){{foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val(cm,0,0,0)*Delta + 0.);end_foreach();}}else if(is_constant(fm.x) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*val(cm,0,0,0)*Delta + 0.);end_foreach();}}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*_const_cm*Delta + 0.);end_foreach();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*_const_cm*Delta + 0.);end_foreach();}}
}





double change (scalar s, scalar sn)
{
  double max = 0.;
  foreach_stencil() {
_stencil_val(cm,0,0,0); {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    }
       
    
#line 302
_stencil_val_a(sn,0,0,0); _stencil_val(s,0,0,0); 
  }end_foreach_stencil();
  
#line 296
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*val(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*_const_cm) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, name))
 return s;}}
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;}}
  }
  return (vector){{-1}};
}
#line 340 "/home/lennard/basilisk/src/utils.h"
#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  if (!(norm > 0.)) qassert ("/home/lennard/basilisk/src/utils.h", 343, "norm > 0.");\
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;\
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -\
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;\
  foreach()\
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {\
      coord _o = {x,y}, _p[2];\
      int _n = 0;\
 if (t.x)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].x = _o.x + _i*Delta/2.;\
     double a = (_p[_n].x - (_S)[0].x)/t.x;\
     _p[_n].y = (_S)[0].y + a*t.y;\
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;\
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&\
    fabs(_p[_n].y - _o.y) <= Delta/2.)\
  _n++;\
     }\
   }\
\
 if (t.y)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].y = _o.y + _i*Delta/2.;\
     double a = (_p[_n].y - (_S)[0].y)/t.y;\
     _p[_n].x = (_S)[0].x + a*t.x;\
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;\
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&\
    fabs(_p[_n].x - _o.x) <= Delta/2.)\
  _n++;\
     }\
   }\
\
      if (_n == 2) {\

#line 380

#line 410 "/home/lennard/basilisk/src/utils.h"
#define end_foreach_segment() } } end_foreach(); }




void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (ferr, " %s", _attribute[s.i].name);}}
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }}}
}

#line 1 "output.h"
#line 1 "/home/lennard/basilisk/src/output.h"
#line 37 "/home/lennard/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};

     
void output_field (struct OutputField p)
{tracing("output_field","/home/lennard/basilisk/src/output.h",46);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  p.n++;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  boundary_internal ((scalar *)p.list, "/home/lennard/basilisk/src/output.h", 58);
  int len = list_len(p.list);
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n - 1);
  int ny = (p.box[1][1] - p.box[0][1])/Delta + 1;
  double ** field = (double **) matrix_new (p.n, ny, len*sizeof(double));
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + p.box[0][0];
    for (int j = 0; j < ny; j++) {
      double y = Delta*j + p.box[0][1];
      if (p.linear) {
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y});}}
      }
      else {
 Point point = locate ((struct _locate){x, y});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : HUGE;}}
      }
    }
  }

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);}}
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + p.box[0][0];
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + p.box[0][1];

 fprintf (p.fp, "%g %g", x, y);
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   fprintf (p.fp, " %g", field[i][len*j + k++]);}}
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (field[0], NULL, len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
end_tracing("output_field","/home/lennard/basilisk/src/output.h",113);}
#line 141 "/home/lennard/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

     
void output_matrix (struct OutputMatrix p)
{tracing("output_matrix","/home/lennard/basilisk/src/output.h",149);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  if (p.linear) {
    scalar f = p.f;
    boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/home/lennard/basilisk/src/output.h", 155);
  }
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 if (!(point.level >= 0)) qassert ("/home/lennard/basilisk/src/output.h", 173, "point.level >= 0");
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
end_tracing("output_matrix","/home/lennard/basilisk/src/output.h",180);}
#line 189 "/home/lennard/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == HUGE) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  if (!(i < 127 - 1)) qassert ("/home/lennard/basilisk/src/output.h", 321, "i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 340 "/home/lennard/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);



    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/home/lennard/basilisk/src/output.h", 422, "pid() == 0");
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  if (!(pid() == 0)) qassert ("/home/lennard/basilisk/src/output.h", 496, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 571 "/home/lennard/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};

     
void output_ppm (struct OutputPPM p)
{tracing("output_ppm","/home/lennard/basilisk/src/output.h",585);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    if (p.spread < 0.)
      p.min = s.min, p.max = s.max;
    else {
      double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
      p.min = avg - spread; p.max = avg + spread;
    }
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)((scalar[]){f, mask,{-1}}), "/home/lennard/basilisk/src/output.h", 608);
    else
      boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/home/lennard/basilisk/src/output.h", 610);
  }

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  if (ny % 2) ny++;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = HUGE;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = HUGE;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
     v = point.level >= 0 ? val(p.f,0,0,0) : HUGE;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = fout;
    if (p.file)
      p.fp = open_image (p.file, p.opt);

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
end_tracing("output_ppm","/home/lennard/basilisk/src/output.h",678);}
#line 710 "/home/lennard/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};

     
void output_grd (struct OutputGRD p)
{tracing("output_grd","/home/lennard/basilisk/src/output.h",720);

  if (!p.fp) p.fp = fout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)((scalar[]){f, mask,{-1}}), "/home/lennard/basilisk/src/output.h", 733);
    else
      boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/home/lennard/basilisk/src/output.h", 735);
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = HUGE;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = HUGE;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 v = point.level >= 0 ? val(p.f,0,0,0) : HUGE;
      }
      if (v == HUGE)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
end_tracing("output_grd","/home/lennard/basilisk/src/output.h",786);}
#line 813 "/home/lennard/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

     
void output_gfs (struct OutputGfs p)
{tracing("output_gfs","/home/lennard/basilisk/src/output.h",842);
  char * fname = p.file;

#if _MPI



  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = fout;
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);




  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if _MPI
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].name)
      cell_size += sizeof(double);}}
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



  {foreach_cell() {
#if _MPI
    if (is_local(cell))
#endif
    {
#if _MPI
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :



      child.x == -1 && child.y == -1 ? 0 :
 child.x == -1 && child.y == 1 ? 1 :
 child.x == 1 && child.y == -1 ? 2 :
 3;
#line 951 "/home/lennard/basilisk/src/output.h"
      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != HUGE ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != HUGE ? - val(s,0,0,0) : (double) DBL_MAX;
     }





   }
   else
     a = is_local(cell) && val(s,0,0,0) != HUGE ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }}}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

#if _MPI
  delete (((scalar[]){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);

#if _MPI
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = fout;
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
end_tracing("output_gfs","/home/lennard/basilisk/src/output.h",1019);}
#line 1043 "/home/lennard/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
  bool unbuffered;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar[]){cm,{-1}}), NULL);
  {scalar*_i=(scalar*)( lista);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);}}
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }}}
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !_MPI
     
void dump (struct Dump p)
{tracing("dump","/home/lennard/basilisk/src/output.h",1096);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,__LINE__);
    strcpy (name, file);
    if (!p.unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  if (!(fp)) qassert ("/home/lennard/basilisk/src/output.h", 1112, "fp");

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar  size=new_scalar("size");
  scalar * list = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

  {foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  pfree (list,__func__,__FILE__,__LINE__);
  if (file) {
    fclose (fp);
    if (!p.unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,__LINE__);
  }delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/lennard/basilisk/src/output.h",1145);}
#else
     
void dump (struct Dump p)
{tracing("dump","/home/lennard/basilisk/src/output.h",1148);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!p.unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar  size=new_scalar("size");
  scalar * list = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);}}
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

  {foreach_cell() {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);}}
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  delete (((scalar[]){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  fclose (fh);
  if (!p.unbuffered && pid() == 0)
    rename (name, file);delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/lennard/basilisk/src/output.h",1219);}
#endif

     
bool restore (struct Dump p)
{tracing("restore","/home/lennard/basilisk/src/output.h",1223);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    {end_tracing("restore","/home/lennard/basilisk/src/output.h",1228);return false;}
  if (!(fp)) qassert ("/home/lennard/basilisk/src/output.h", 1229, "fp");

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }
#line 1260 "/home/lennard/basilisk/src/output.h"
  init_grid (1 << header.depth);



  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }}}
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,__LINE__);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }
#line 1339 "/home/lennard/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1}}});



  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }

    fseek (fp, sizeof(double), SEEK_CUR);
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX)
 val(s,0,0,0) = val;
    }}}
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}


  scalar * other = NULL;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);}}
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  {end_tracing("restore","/home/lennard/basilisk/src/output.h",1389);return true;}
end_tracing("restore","/home/lennard/basilisk/src/output.h",1390);}
#line 431 "/home/lennard/basilisk/src/utils.h"
#line 4 "/home/lennard/basilisk/src/predictor-corrector.h"



extern scalar * evolving;

double (* update) (scalar * evolving, scalar * updates, double dtmax) = NULL;



double (* gradient) (double, double, double) = minmod2;


double dt = 0.;

     
static void advance_generic (scalar * output, scalar * input, scalar * updates,
        double dt)
{tracing("advance_generic","/home/lennard/basilisk/src/predictor-corrector.h",19);
  if (input != output)
    trash (output);
  foreach_stencil() {
    scalar o, i, u;
    {scalar*_i0=updates;scalar*_i1=input;scalar*_i2= output;if(_i0)for(u=*_i0,i=*_i1,o=*_i2;_i0->i>= 0;u=*++_i0,i=*++_i1,o=*++_i2){
      {_stencil_val_a(o,0,0,0); _stencil_val(i,0,0,0);_stencil_val(u,0,0,0);   }}}
  }end_foreach_stencil();
  {
#line 24
foreach() {
    scalar o, i, u;
    {scalar*_i0=updates;scalar*_i1=input;scalar*_i2= output;if(_i0)for(u=*_i0,i=*_i1,o=*_i2;_i0->i>= 0;u=*++_i0,i=*++_i1,o=*++_i2){
      val(o,0,0,0) = val(i,0,0,0) + dt*val(u,0,0,0);}}
  }end_foreach();}
end_tracing("advance_generic","/home/lennard/basilisk/src/predictor-corrector.h",29);}

static void (* advance) (scalar * output, scalar * input, scalar * updates,
    double dt) = advance_generic;

static int defaults_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults(const int i,const double t,Event *_ev){tracing("defaults","/home/lennard/basilisk/src/predictor-corrector.h",34);
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].gradient = gradient;}}

  display ((struct _display){"box();"});
}{end_tracing("defaults","/home/lennard/basilisk/src/predictor-corrector.h",41);return 0;}end_tracing("defaults","/home/lennard/basilisk/src/predictor-corrector.h",41);}

     
void run()
{tracing("run","/home/lennard/basilisk/src/predictor-corrector.h",44);
  t = 0., iter = 0;
  init_grid (N);


  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {

    scalar * updates = list_clone (evolving);
    dt = dtnext (update (evolving, updates, DT));
    if (gradient != zero) {

      scalar * predictor = list_clone (evolving);

      advance (predictor, evolving, updates, dt/2.);

      update (predictor, updates, dt);
      delete (predictor);
      pfree (predictor,__func__,__FILE__,__LINE__);
    }
    advance (evolving, evolving, updates, dt);
    delete (updates);
    pfree (updates,__func__,__FILE__,__LINE__);
    update_perf();
    iter = inext, t = tnext;
  }
  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
end_tracing("run","/home/lennard/basilisk/src/predictor-corrector.h",75);}
#line 78 "./qg.h"
#line 1 "nodal-poisson.h"
#line 1 "./nodal-poisson.h"
#line 11 "./nodal-poisson.h"
#line 1 "inner-vertex.h"
#line 1 "./inner-vertex.h"



#define foreach_inner_face_generic()\
  OMP_PARALLEL() {\
\
  int index_left = 1, index_right = 1;\
  int index_top = 1, index_bottom = 1;\
  int index_front = 1, index_back = 1;\
\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2 + index_left ; _k <= point.n + 2 - index_right; _k++) {\
    point.i = _k;\
\
    for (point.j = 2 + index_bottom; point.j <= point.n + 2 - index_top; point.j++)\
\
\
\
        {\
\
   POINT_VARIABLES\

#line 36

#define end_foreach_inner_face_generic()\
\
 }\
\
  }\
}\

#line 43


#define foreach_inner_vertex()\
foreach_inner_face_generic() {\
  x -= Delta/2.;\
\
  y -= Delta/2.;\
\
\
\
\

#line 54

#define end_foreach_inner_vertex() } end_foreach_inner_face_generic()



#define foreach_inner_vertex_level(l)\
OMP_PARALLEL() {\
\
  int index_left = 1, index_right = 1;\
  int index_top = 1, index_bottom = 1;\
  int index_front = 1, index_back = 1;\
\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2 + index_left; _k <= point.n + 2 - index_right; _k++) {\
    point.i = _k;\
\
    for (point.j = 2 + index_bottom; point.j <= point.n + 2 - index_top; point.j++)\
\
\
\
 {\
\
          POINT_VARIABLES\
  x -= Delta/2.;\
\
  y -= Delta/2.;\
\
\
\
\

#line 98

#define end_foreach_inner_vertex_level()\
\
 }\
\
  }\
}\

#line 105


#define foreach_vertex_level(l)\
OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2 + 1; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2 + 1; point.j++)\
\
\
\
 {\
\
          POINT_VARIABLES\
  x -= Delta/2.;\
\
  y -= Delta/2.;\
\
\
\
\

#line 131

#define end_foreach_vertex_level()\
\
 }\
\
  }\
}\

#line 138

#line 12 "./nodal-poisson.h"
#line 1 "my_vertex.h"
#line 1 "./my_vertex.h"
#line 17 "./my_vertex.h"
#define foreach_vert()\
  update_cache();\
  foreach_cache(tree->leaves) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
\
\

#line 27

#define end_foreach_vert() } end_foreach_cache()

#define foreach_vert_level(l) {\
  if (l <= depth()) {\
    update_cache();\
    CacheLevel _active = tree->active[l];\
    foreach_cache_level (_active,l) {\
      x -= Delta/2.;\
\
      y -= Delta/2.;\
\
\
\
\

#line 42

#define end_foreach_vert_level() } end_foreach_cache_level(); }}





static inline void restriction_vert (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  val(s,0,0,0) = fine(s,0,0,0);
}



static inline void restriction_coarsen_vert (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  val(s,0,0,0) = (fine(s,1,0,0) + 2*fine(s,0,0,0) + fine(s,-1,0,0) +
  fine(s,0,1,0) + fine(s,0,-1,0))/6.;

}





static inline void refine_vert (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  fine (s,0,0,0) = val(s,0,0,0);

  fine (s,1,0,0) = (val(s,0,0,0) + val(s,1,0,0))/2.;

  fine (s,0,1,0) = (val(s,0,0,0) + val(s,0,1,0))/2.;






  fine(s,1,1,0) = (val(s,0,0,0) + val(s,1,0,0) + val(s,0,1,0) + val(s,1,1,0))/4.;
#line 91 "./my_vertex.h"
}


static inline void prolongation_vert (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  refine_vert(point, s);

      fine(s,0,0,0) = (val(s,-1,0,0) + val(s,1,0,0)

         + val(s,0,1,0) + val(s,0,-1,0)

         )/(2*2);
}




static inline void refine_vert4 (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  fine (s,0,0,0) = val(s,0,0,0);

  fine (s,1,0,0) = (9*(val(s,0,0,0) + val(s,1,0,0)) - (val(s,-1,0,0) + val(s,2,0,0)))/16.;

  fine (s,0,1,0) = (9*(val(s,0,0,0) + val(s,0,1,0)) - (val(s,0,-1,0) + val(s,0,2,0)))/16.;






  fine(s,1,1,0) = (9.*(val(s,0,0,0) + val(s,1,1,0) + val(s,0,1,0) + val(s,1,0,0)) -
     (val(s,-1,-1,0) + val(s,2,2,0) + val(s,-1,2,0) + val(s,2,-1,0)))/32.;
#line 133 "./my_vertex.h"
}

static inline void prolongation_vert4 (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  refine_vert4(point, s);
  if (!is_leaf(cell)) {
    fine(s,0,0,0) = (4*(fine(s,-1,0,0) + fine(s,1,0,0)) -
       (fine(s,-2,0,0) + fine(s,2,0,0))

       +4* (fine(s,0,-1,0) + fine(s,0,1,0)) -
       (fine(s,0,-2,0) + fine(s,0,2,0))

       )/(2*6);
  }
}


double interp4_x (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return (4*(val(s,-1,0,0) + val(s,1,0,0))- (val(s,-2,0,0) + val(s,2,0,0)))/6.;
}

#line 149
double interp4_y (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return (4*(val(s,0,-1,0) + val(s,0,1,0))- (val(s,0,-2,0) + val(s,0,2,0)))/6.;
}




double interp5_x (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return ((3./128)*coarse(s,-2,0,0) - (5./32)*coarse(s,-1,0,0) +
   (45./64.)*coarse(s,0,0,0) + (15./32.)*coarse(s,1,0,0) -
   (5./128.)*coarse(s,2,0,0));
}

double interp5_y (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return ((3./128)*coarse(s,0,-2,0) - (5./32)*coarse(s,0,-1,0) +
   (45./64.)*coarse(s,0,0,0) + (15./32.)*coarse(s,0,1,0) -
   (5./128.)*coarse(s,0,2,0));
}

static inline void refine_vert5 (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  fine (s,0,0,0) = val(s,0,0,0);

  fine (s,1,0,0) = ((3./128)*val(s,-2,0,0) - (5./32)*val(s,-1,0,0) + (45./64.)*val(s,0,0,0)
      + (15./32.)*val(s,1,0,0) - (5./128.)*val(s,2,0,0));

  fine (s,0,1,0) = ((3./128)*val(s,0,-2,0) - (5./32)*val(s,0,-1,0) + (45./64.)*val(s,0,0,0)
      + (15./32.)*val(s,0,1,0) - (5./128.)*val(s,0,2,0));







  fine (s,1,1,0) = ((3./128)*val(s,-2,-2,0) - (5./32)*val(s,-1,-1,0) + (45./64.)*val(s,0,0,0) +
      (15./32.)*val(s,1,1,0) - (5./128.)*val(s,2,2,0));
#line 196 "./my_vertex.h"
}
#line 13 "./nodal-poisson.h"
#line 1 "poisson.h"
#line 1 "/home/lennard/basilisk/src/poisson.h"
#line 32 "/home/lennard/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
      {foreach_level_or_leaf (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = 0.;}}end_foreach_level_or_leaf();}





    else
      {foreach_level (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = bilinear (point, s);}}end_foreach_level();}





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




  foreach_stencil() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 {_stencil_val_r(s,0,0,0); _stencil_val(ds,0,0,0); }}}
  }end_foreach_stencil();




  {
#line 84
foreach() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 val(s,0,0,0) += val(ds,0,0,0);}}
  }end_foreach();}
}
#line 102 "/home/lennard/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 125 "/home/lennard/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);






  for (int b = 0; b < nboundary; b++)
    {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];}}




  mgstats s = {0};
  double sum = 0.;
  foreach_stencil ()
    {scalar*_i=(scalar*)( p.b);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      { _stencil_val(s,0,0,0); }}}end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)){
#line 164
foreach ()
    {scalar*_i=(scalar*)( p.b);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      sum += val(s,0,0,0);}}end_foreach();mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 167
s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data,
       s.nrelax,
       p.minlevel,
       grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);
#line 199 "/home/lennard/basilisk/src/poisson.h"
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = p.minlevel;




  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 258 "/home/lennard/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
          vector alpha;
          scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;



};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
#line 296 "/home/lennard/basilisk/src/poisson.h"
  scalar c = a;






  if(!is_constant(lambda) && !is_constant(alpha.x)){{foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 305
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }
#line 319 "/home/lennard/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 305
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }
#line 319 "/home/lennard/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(!is_constant(lambda) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 305
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }
#line 319 "/home/lennard/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 305
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }
#line 319 "/home/lennard/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}
#line 338 "/home/lennard/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
  double maxres = 0.;
#line 372 "/home/lennard/basilisk/src/poisson.h"
  foreach_stencil () {
    _stencil_val_a(res,0,0,0); _stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);  
    
      {_stencil_val_r(res,0,0,0);_stencil_val(alpha.x,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
  _stencil_val(alpha.x,1,0,0);_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);     }
      
#line 375
{_stencil_val_r(res,0,0,0);_stencil_val(alpha.y,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
  _stencil_val(alpha.y,0,1,0);_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);     }






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 385
}end_foreach_stencil();
#line 372 "/home/lennard/basilisk/src/poisson.h"
  if(!is_constant(lambda) && !is_constant(alpha.x)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 372
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) += (val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 375
val(res,0,0,0) += (val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 385
}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);
#line 372 "/home/lennard/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 372
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) += (val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 375
val(res,0,0,0) += (val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 385
}else if(!is_constant(lambda) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
#line 372 "/home/lennard/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 372
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) += (_const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  _const_alpha.x*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 375
val(res,0,0,0) += (_const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  _const_alpha.y*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 385
}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
#line 372 "/home/lennard/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 372
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) += (_const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  _const_alpha.x*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 375
val(res,0,0,0) += (_const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  _const_alpha.y*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 385
}

  return maxres;
}
#line 399 "/home/lennard/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar[]){alpha.x,alpha.y,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;




  mgstats s = mg_solve ((struct MGSolve){((scalar[]){a,{-1}}), ((scalar[]){b,{-1}}), residual, relax,
   &p, p.nrelax, p.res, .minlevel = max(1, p.minlevel)});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 461 "/home/lennard/basilisk/src/poisson.h"
struct Project {
  vector uf;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};

     
mgstats project (struct Project q)
{tracing("project","/home/lennard/basilisk/src/poisson.h",470);
  vector uf = q.uf;
  scalar p = q.p;
          vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar  div=new_scalar("div");
  foreach_stencil() {
    _stencil_val_a(div,0,0,0);  
    
      {_stencil_val_r(div,0,0,0); _stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);  }
      
#line 487
{_stencil_val_r(div,0,0,0); _stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);  }
    _stencil_val_r(div,0,0,0);  
  }end_foreach_stencil();
  {
#line 484
foreach() {
    val(div,0,0,0) = 0.;
    
      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
      
#line 487
val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);
    val(div,0,0,0) /= dt*Delta;
  }end_foreach();}
#line 500 "/home/lennard/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_r(uf.x,0,0,0);_stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);   }}end__stencil_is_face_x()
#line 506
_stencil_is_face_y(){
    {_stencil_val_r(uf.y,0,0,0);_stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();




  
#line 506
if(!is_constant(alpha.x)){{foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 506
is_face_y(){
    val(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);




  {
#line 506
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*_const_alpha.x*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 506
is_face_y(){
    val(uf.y,0,0,0) -= dt*_const_alpha.y*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}

  {delete((scalar*)((scalar[]){div,{-1}}));{end_tracing("project","/home/lennard/basilisk/src/poisson.h",509);return mgp;}}delete((scalar*)((scalar[]){div,{-1}}));
end_tracing("project","/home/lennard/basilisk/src/poisson.h",510);}static double _boundary0(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
#line 26 "./nodal-poisson.h"
0.;}}static double _boundary1(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0.;}}static double _boundary2(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0.;}}static double _boundary3(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0.;}}static double _boundary4(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
#line 56
0.;}}static double _boundary5(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0.;}}static double _boundary6(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0.;}}static double _boundary7(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0.;}}
#line 14 "./nodal-poisson.h"

mgstats vpoisson (struct Poisson p) {

  scalar  da=new_vertex_scalar("da"),  res=new_vertex_scalar("res"), a = p.a, b = p.b;
  scalar_clone (da, a);
  _attribute[da.i].restriction = restriction_vert;
  _attribute[da.i].prolongation = refine_vert;

  scalar  mask=new_vertex_scalar("mask");
  _attribute[mask.i].restriction = restriction_vert;
  _attribute[mask.i].prolongation = refine_vert;

_attribute[mask.i].dirty=1,_attribute[mask.i].boundary[left]=_boundary0,_attribute[mask.i].boundary_homogeneous[left]=_boundary0;
_attribute[mask.i].dirty=1,_attribute[mask.i].boundary[right]=_boundary1,_attribute[mask.i].boundary_homogeneous[right]=_boundary1;
_attribute[mask.i].dirty=1,_attribute[mask.i].boundary[top]=_boundary2,_attribute[mask.i].boundary_homogeneous[top]=_boundary2;
_attribute[mask.i].dirty=1,_attribute[mask.i].boundary[bottom]=_boundary3,_attribute[mask.i].boundary_homogeneous[bottom]=_boundary3;







  foreach_vertex_stencil()
    {_stencil_val_a(mask,0,0,0);  }end_foreach_vertex_stencil();







  {
#line 37
foreach_vertex()
    val(mask,0,0,0) = 1.;end_foreach_vertex();}
  boundary_internal ((scalar *)((scalar[]){mask,{-1}}), "./nodal-poisson.h", 39);
  restriction(((scalar[]){mask,{-1}}));
  for (int l = 0; l <= depth(); l++) {
    boundary_level(((scalar[]){mask,{-1}}), l);

  }

  if (p.res)
    res = p.res[0];
  else
    scalar_clone (res, b);

  _attribute[res.i].prolongation = refine_vert;




_attribute[res.i].dirty=1,_attribute[res.i].boundary[left]=_boundary4,_attribute[res.i].boundary_homogeneous[left]=_boundary4;
_attribute[res.i].dirty=1,_attribute[res.i].boundary[right]=_boundary5,_attribute[res.i].boundary_homogeneous[right]=_boundary5;
_attribute[res.i].dirty=1,_attribute[res.i].boundary[top]=_boundary6,_attribute[res.i].boundary_homogeneous[top]=_boundary6;
_attribute[res.i].dirty=1,_attribute[res.i].boundary[bottom]=_boundary7,_attribute[res.i].boundary_homogeneous[bottom]=_boundary7;


  mgstats mg; mg.sum = HUGE, mg.resa = HUGE;
  mg.nrelax = p.nrelax ? p.nrelax : 5;
  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;
  mg.minlevel = p.minlevel ? p.minlevel : 1;

  for (mg.i = 0; mg.i < NITERMAX; mg.i++) {

    double max = 0;
    foreach_vertex_stencil() {
      _stencil_val_a(res,0,0,0); _stencil_val(b,0,0,0);_stencil_val(mask,0,0,0); 
       {
          _stencil_val_r(res,0,0,0);_stencil_val(a,-1,0,0);_stencil_val(a,0,0,0); _stencil_val(a,1,0,0);_stencil_val(mask,0,0,0);     
      } 
#line 74
{
          _stencil_val_r(res,0,0,0);_stencil_val(a,0,-1,0);_stencil_val(a,0,0,0); _stencil_val(a,0,1,0);_stencil_val(mask,0,0,0);     
      }
_stencil_val(res,0,0,0);
        {_stencil_val(res,0,0,0);  }
         
    
#line 79
}end_foreach_vertex_stencil();
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction (max:max)){
#line 72
foreach_vertex() {
      val(res,0,0,0) = val(b,0,0,0)*val(mask,0,0,0);
       {
          val(res,0,0,0) -= (val(a,-1,0,0) - 2.*val(a,0,0,0) + val(a,1,0,0))/(sq(Delta))*val(mask,0,0,0);
      } 
#line 74
{
          val(res,0,0,0) -= (val(a,0,-1,0) - 2.*val(a,0,0,0) + val(a,0,1,0))/(sq(Delta))*val(mask,0,0,0);
      }
      if (fabs(val(res,0,0,0)) > max)
        max = fabs(val(res,0,0,0));
    }end_foreach_vertex();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

    
#line 81
mg.resa = max;
    if (mg.i == 0)
      mg.resb = max;

    if (max < TOLERANCE && mg.i >= NITERMIN)
      break;

    _attribute[res.i].restriction = restriction_vert;
    boundary_internal ((scalar *)((scalar[]){res,{-1}}), "./nodal-poisson.h", 89);
    _attribute[res.i].restriction = restriction_coarsen_vert;
    multigrid_restriction (((scalar[]){res,{-1}}));
    for (int l = 0; l <= depth(); l++) {
      boundary_level(((scalar[]){res,{-1}}), l);
    }


    {foreach_vertex_level(mg.minlevel)
      val(da,0,0,0) = 0;end_foreach_vertex_level();}

    for (int l = mg.minlevel; l <= depth(); l++) {
      boundary_level(((scalar[]){da,{-1}}), l);

      for (int rel = 0; rel < mg.nrelax; rel++) {
        {foreach_vertex_level(l) {
          double d = 0;
          val(da,0,0,0) = -val(res,0,0,0)*sq(Delta);
           {
              val(da,0,0,0) += (val(da,1,0,0) + val(da,-1,0,0))*val(mask,0,0,0);
              d += 2.;
          } 
#line 107
{
              val(da,0,0,0) += (val(da,0,1,0) + val(da,0,-1,0))*val(mask,0,0,0);
              d += 2.;
          }
          val(da,0,0,0) /= d;
        }end_foreach_vertex_level();}
        boundary_level(((scalar[]){da,{-1}}), l);
      }





      if (l < depth()){
        {foreach_vertex_level(l){

          refine_vert (point, da);
        }end_foreach_vertex_level();}
        boundary_level(((scalar[]){da,{-1}}), l+1);

      }
    }

    foreach_vertex_stencil()
      {_stencil_val_r(a,0,0,0); _stencil_val(da,0,0,0); }end_foreach_vertex_stencil();

    {
#line 130
foreach_vertex()
      val(a,0,0,0) += val(da,0,0,0);end_foreach_vertex();}
    boundary_internal ((scalar *)((scalar[]){a,{-1}}), "./nodal-poisson.h", 132);
  }
  if (mg.resa > TOLERANCE)
    fprintf (ferr, "Convergence for %s not reached.\n"
      "mg.i = %d, mg.resb: %g mg.resa: %g\n",
      _attribute[a.i].name, mg.i, mg.resb, mg.resa);
  if (p.tolerance)
    TOLERANCE = defaultol;
  {delete((scalar*)((scalar[]){mask,res,da,{-1}}));return mg;}delete((scalar*)((scalar[]){mask,res,da,{-1}}));
}
#line 79 "./qg.h"





double Re = 0.;
double delta_nl = 0.;

double iend = 0.;

double f0 = 1.0;
double beta = 0.;
double hEkb = 0.;
double tau0 = 0.;
double nu = 0.;
double sbc = 0.;
double tend = 100.;
double dtout = 1.;
double dh[1] = {1.};

char dpath[80];




scalar  psi={0};
scalar  q={1};
scalar  zeta={2};





double bc_fac = 0;
scalar * evolving = NULL;
mgstats mgpsi;static double _boundary8(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
#line 153
0;}}static double _boundary9(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0;}}static double _boundary10(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0;}}static double _boundary11(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
0;}}static double _boundary12(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 

bc_fac/sq(Delta)*(val(psi,1,0,0) - val(psi,0,0,0));}}static double _boundary13(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
bc_fac/sq(Delta)*(val(psi,0,0,0) - val(psi,1,0,0));}}static double _boundary14(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
bc_fac/sq(Delta)*(val(psi,0,1,0) - val(psi,0,0,0));}}static double _boundary15(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
bc_fac/sq(Delta)*(val(psi,0,0,0) - val(psi,0,1,0));}}static double _boundary16(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 

bc_fac/sq(Delta)*(val(psi,1,0,0) - val(psi,0,0,0));}}static double _boundary17(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
bc_fac/sq(Delta)*(val(psi,0,0,0) - val(psi,1,0,0));}}static double _boundary18(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
bc_fac/sq(Delta)*(val(psi,0,1,0) - val(psi,0,0,0));}}static double _boundary19(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return 
bc_fac/sq(Delta)*(val(psi,0,0,0) - val(psi,0,1,0));}}
#line 148 "./qg.h"
void set_bc()
{

  bc_fac = sbc/((0.5*sbc + 1));

_attribute[psi.i].dirty=1,_attribute[psi.i].boundary[left]=_boundary8,_attribute[psi.i].boundary_homogeneous[left]=_boundary8;
_attribute[psi.i].dirty=1,_attribute[psi.i].boundary[right]=_boundary9,_attribute[psi.i].boundary_homogeneous[right]=_boundary9;
_attribute[psi.i].dirty=1,_attribute[psi.i].boundary[bottom]=_boundary10,_attribute[psi.i].boundary_homogeneous[bottom]=_boundary10;
_attribute[psi.i].dirty=1,_attribute[psi.i].boundary[top]=_boundary11,_attribute[psi.i].boundary_homogeneous[top]=_boundary11;

_attribute[zeta.i].dirty=1,_attribute[zeta.i].boundary[left]=_boundary12,_attribute[zeta.i].boundary_homogeneous[left]=_boundary12;
_attribute[zeta.i].dirty=1,_attribute[zeta.i].boundary[right]=_boundary13,_attribute[zeta.i].boundary_homogeneous[right]=_boundary13;
_attribute[zeta.i].dirty=1,_attribute[zeta.i].boundary[bottom]=_boundary14,_attribute[zeta.i].boundary_homogeneous[bottom]=_boundary14;
_attribute[zeta.i].dirty=1,_attribute[zeta.i].boundary[top]=_boundary15,_attribute[zeta.i].boundary_homogeneous[top]=_boundary15;

_attribute[q.i].dirty=1,_attribute[q.i].boundary[left]=_boundary16,_attribute[q.i].boundary_homogeneous[left]=_boundary16;
_attribute[q.i].dirty=1,_attribute[q.i].boundary[right]=_boundary17,_attribute[q.i].boundary_homogeneous[right]=_boundary17;
_attribute[q.i].dirty=1,_attribute[q.i].boundary[bottom]=_boundary18,_attribute[q.i].boundary_homogeneous[bottom]=_boundary18;
_attribute[q.i].dirty=1,_attribute[q.i].boundary[top]=_boundary19,_attribute[q.i].boundary_homogeneous[top]=_boundary19;

}




     
void invertq(scalar psi, scalar q)
{tracing("invertq","./qg.h",174);

  mgpsi = vpoisson((struct Poisson){psi, q});
  boundary_internal ((scalar *)((scalar[]){psi,{-1}}), "./qg.h", 178);

  set_bc();
  boundary_internal ((scalar *)((scalar[]){q,{-1}}), "./qg.h", 181);
end_tracing("invertq","./qg.h",182);}

     
void comp_del2(scalar psi, scalar zeta, double add, double fac)
{tracing("comp_del2","./qg.h",185);
  {foreach_inner_vertex()
      val(zeta,0,0,0) = add*val(zeta,0,0,0) + fac*(val(psi,1,0,0) + val(psi,-1,0,0) + val(psi,0,1,0) + val(psi,0,-1,0) - 4*val(psi,0,0,0))/(sq(Delta));end_foreach_inner_vertex();}

  boundary_internal ((scalar *)((scalar[]){zeta,{-1}}), "./qg.h", 190);
end_tracing("comp_del2","./qg.h",191);}

     
void comp_q(scalar psi, scalar q)
{tracing("comp_q","./qg.h",194);
  comp_del2 (psi, q, 0., 1.);

  boundary_internal ((scalar *)((scalar[]){q,{-1}}), "./qg.h", 198);
end_tracing("comp_q","./qg.h",199);}




     
double advection_pv(scalar zeta, scalar q, scalar psi, scalar dqdt, double dtmax)
{tracing("advection_pv","./qg.h",205);
  {foreach_inner_vertex()
    val(dqdt,0,0,0) += -((( val(psi,1, 0,0 )-val(psi,-1, 0,0))*(val(zeta,0, 1,0)-val(zeta, 0 ,-1,0)) +(val(psi,0 ,-1,0)-val(psi, 0 ,1,0))*(val(zeta,1, 0,0)-val(zeta,-1, 0,0 )) + val(psi,1, 0,0 )*( val(zeta,1,1,0 ) - val(zeta,1,-1,0 )) - val(psi,-1, 0,0)*( val(zeta,-1,1,0) - val(zeta,-1,-1,0)) - val(psi, 0 ,1,0)*( val(zeta,1,1,0 ) - val(zeta,-1,1,0 )) + val(psi,0 ,-1,0)*( val(zeta,1,-1,0) - val(zeta,-1,-1,0)) + val(zeta, 0 ,1,0)*( val(psi,1,1,0 ) - val(psi,-1,1,0 )) - val(zeta,0 ,-1,0)*( val(psi,1,-1,0) - val(psi,-1,-1,0)) - val(zeta,1, 0,0 )*( val(psi,1,1,0 ) - val(psi,1,-1,0 )) + val(zeta,-1, 0,0)*( val(psi,-1,1,0) - val(psi,-1,-1,0))) /(12.*Delta*Delta)) - (beta*(val(psi,1,0,0) - val(psi,-1,0,0))/(2*Delta));end_foreach_inner_vertex();}


  static double previous = 0.;
  dtmax /= CFL;
  foreach_face_stencil(){_stencil_is_face_x(){{    
       _stencil_val(psi,0,0,0);_stencil_val(psi,0,1,0);   
           
           
       
         
  }}end__stencil_is_face_x()
#line 213
_stencil_is_face_y(){{    
       _stencil_val(psi,0,0,0);_stencil_val(psi,1,0,0);   
           
           
       
         
  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 213
foreach_face_generic(){is_face_x(){{
      double u = (val(psi,0,1,0) - val(psi,0,0,0))/Delta;
      if (u != 0.) {
        double dt = Delta/fabs(u);
        if (dt < dtmax) dtmax = dt;
      }
  }}end_is_face_x()
#line 213
is_face_y(){{
      double u = (val(psi,1,0,0) - val(psi,0,0,0))/Delta;
      if (u != 0.) {
        double dt = Delta/fabs(u);
        if (dt < dtmax) dtmax = dt;
      }
  }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

  
#line 221
dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  {end_tracing("advection_pv","./qg.h",225);return dtmax;}
end_tracing("advection_pv","./qg.h",226);}

     
void dissip (scalar zeta, scalar dqdt)
{tracing("dissip","./qg.h",229);
  comp_del2(zeta, dqdt, 1., nu);
end_tracing("dissip","./qg.h",232);}





     
void ekman_friction (scalar zeta, scalar dqdt)
{tracing("ekman_friction","./qg.h",239);
  {foreach_inner_vertex()
    val(dqdt,0,0,0) -= hEkb*f0/(2*dh[0])*val(zeta,0,0,0);end_foreach_inner_vertex();}
end_tracing("ekman_friction","./qg.h",243);}




     
void surface_forcing (scalar dqdt)
{tracing("surface_forcing","./qg.h",249);
  {foreach_inner_vertex()
    val(dqdt,0,0,0) -= tau0/L0*3.14159265358979*sin(3.14159265358979*y/L0);end_foreach_inner_vertex();}
end_tracing("surface_forcing","./qg.h",253);}


static void advance_qg (scalar * output, scalar * input,
                        scalar * updates, double dt)
{

  scalar qi = input[0];
  scalar qo = output[0];
  scalar dq = updates[0];

  {foreach_inner_vertex()
    val(qo,0,0,0) = val(qi,0,0,0) + val(dq,0,0,0)*dt;end_foreach_inner_vertex();}

}


double update_qg (scalar * evolving, scalar * updates, double dtmax)
{


  foreach_vertex_stencil()
    {scalar*_i=(scalar*)( updates);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
        {_stencil_val_a(s,0,0,0);  }}}end_foreach_vertex_stencil();


  {
#line 274
foreach_vertex()
    {scalar*_i=(scalar*)( updates);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
        val(s,0,0,0) = 0.;}}end_foreach_vertex();}

  scalar q = evolving[0];
  scalar dqdt = updates[0];

  invertq(psi, q);
  comp_del2(psi, zeta, 0., 1.0);
  dtmax = advection_pv(zeta, q, psi, dqdt, dtmax);

  dissip(zeta, dqdt);
  ekman_friction(zeta, dqdt);
  surface_forcing(dqdt);

  return dtmax;
}

void set_vars()
{

  _attribute[psi.i].restriction = restriction_vert;
  _attribute[psi.i].prolongation = refine_vert;

  _attribute[zeta.i].restriction = restriction_vert;
  _attribute[zeta.i].prolongation = refine_vert;

  _attribute[q.i].restriction = restriction_vert;
  _attribute[q.i].prolongation = refine_vert;


  reset (((scalar[]){q, psi, zeta,{-1}}), 0.);
#line 316 "./qg.h"
  evolving = list_copy(((scalar[]){q,{-1}}));
  advance = advance_qg;
  update = update_qg;
  fprintf(fout,"ok\n");
}


static int defaults_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults_0(const int i,const double t,Event *_ev){tracing("defaults_0","./qg.h",323);{
  set_bc();
  set_vars();
}{end_tracing("defaults_0","./qg.h",326);return 0;}end_tracing("defaults_0","./qg.h",326);}





void set_const() {
  comp_q(psi,q);
  boundary_internal ((scalar *)all, "./qg.h", 334);
}

static int init_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int init(const int i,const double t,Event *_ev){tracing("init","./qg.h",337); {
  set_const();
}{end_tracing("init","./qg.h",339);return 0;}end_tracing("init","./qg.h",339);}




void trash_vars(){

  pfree(evolving,__func__,__FILE__,__LINE__);


}

static int cleanup_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 1234567890);*ip=i;*tp=t;return ret;}      static int cleanup(const int i,const double t,Event *_ev){tracing("cleanup","./qg.h",351);
{
  trash_vars();
}{end_tracing("cleanup","./qg.h",354);return 0;}end_tracing("cleanup","./qg.h",354);}
#line 22 "qg.c"
#line 1 "extra.h"
#line 1 "./extra.h"







#line 1 "/home/lennard/basilisk/src/ast/std/sys/stat.h"
#include <sys/stat.h>
#line 9 "./extra.h"
#line 1 "/home/lennard/basilisk/src/ast/std/sys/types.h"
#include <sys/types.h>
#line 10 "./extra.h"

void trim_whitespace(char* s) {
  const char* d = s;
  do {
    while (*d == ' ')
      ++d;
  } while (*s++ = *d++);
  *s = '\0';
}


void str2array(char *tmps2, double *array){
  char* p;
  int n = 0;
  p = strtok(tmps2,"[,]");
  while (p != NULL){
    array[n] = atof(p);
    p = strtok(NULL, ",");
    n += 1;
  }
}

void read_params(char* path2file)
{
  FILE * fp;
  if ((fp = fopen(path2file, "rt"))) {
    char tempbuff[300];
    while(fgets(tempbuff,300,fp)) {
      trim_whitespace(tempbuff);
      char* tmps1 = strtok(tempbuff, "=");
      char* tmps2 = strtok(NULL, "=");

      if (strcmp(tmps1,"N") ==0) { N = atoi(tmps2); }

      else if (strcmp(tmps1,"DT") ==0) { DT = atof(tmps2); }
      else if (strcmp(tmps1,"CFL") ==0) { CFL = atof(tmps2); }
      else if (strcmp(tmps1,"TOLERANCE")==0) { TOLERANCE= atof(tmps2); }

      else if (strcmp(tmps1,"delta_nl")==0){
        delta_nl = atof(tmps2);
        L0 = pow(tau0/(pow(delta_nl*beta,2.)), 1./3.);
        fprintf(fout, "nonlinear thickness is %g.\n", pow(tau0/(pow(L0, 3.)*pow(beta, 2.)),0.5));
      }
      else if (strcmp(tmps1,"Re")==0){
        Re = atof(tmps2);
        nu = pow(tau0*pow(L0,3.)/pow(Re,4.), 1./2.);
        fprintf(fout, "Sublayer Re is %g.\n", pow(tau0*pow(L0,3.)/pow(nu,2.),1./4.));
      }
      else if (strcmp(tmps1,"f0") ==0) { f0 = atof(tmps2); }
      else if (strcmp(tmps1,"beta") ==0) { beta = atof(tmps2); }
      else if (strcmp(tmps1,"hEkb") ==0) { hEkb = atof(tmps2); }
      else if (strcmp(tmps1,"tau0") ==0) { tau0 = atof(tmps2); }
      else if (strcmp(tmps1,"sbc") ==0) { sbc = atof(tmps2); }
      else if (strcmp(tmps1,"tend") ==0) { tend = atof(tmps2); }
      else if (strcmp(tmps1,"dtout")==0) { dtout = atof(tmps2); }
      else if (strcmp(tmps1,"dh") ==0) { str2array(tmps2, dh);}
      else if (strcmp(tmps1,"iend") ==0) { iend = atof(tmps2); }


    }
    fclose(fp);
  } else {
    fprintf(fout, "file %s not found\n", path2file);
    exit(0);
  }




  if (beta != 0) DT = min(3.1415/(2.*beta*L0),0.5*min(DT,sq(L0/N)/nu/4.));
  else 0.5*min(DT,sq(L0/N)/nu/4.);

  fprintf(fout, "Config: N = %d, L0 = %g\n", N, L0);

}




void create_outdir()
{
  if (pid() == 0) {
    for (int i=1; i<10000; i++) {
      sprintf(dpath, "outdir_%04d/", i);
      if (mkdir(dpath, 0777) == 0) {
        fprintf(fout,"Writing output in %s\n",dpath);
        break;
      }
    }
  }
#if _MPI
  MPI_Bcast(&dpath, 80, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
}

void backup_config()
{
  fprintf(fout, "Backup config\n");
  char ch;
  char name[120];
  sprintf (name,"%sparams.in", dpath);
  FILE * source = fopen("params.in", "r");
  FILE * target = fopen(name, "w");
  while ((ch = fgetc(source)) != EOF)
    fputc(ch, target);
  fclose(source);
  fclose(target);
}
#line 23 "qg.c"
#line 1 "netcdf_vertex_bas.h"
#line 1 "./netcdf_vertex_bas.h"




#line 1 "/home/lennard/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 6 "./netcdf_vertex_bas.h"
#line 1 "/home/lennard/basilisk/src/ast/std/string.h"
#include <string.h>
#line 7 "./netcdf_vertex_bas.h"
#line 1 "/usr/include/netcdf.h"
#line 16 "/usr/include/netcdf.h"
#line 1 "/home/lennard/basilisk/src/ast/std/stddef.h"

#line 1 "/home/lennard/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 17 "/usr/include/netcdf.h"
#line 1 "/usr/include/errno.h"
#line 25 "/usr/include/errno.h"
#line 1 "/usr/include/features.h"
#line 392 "/usr/include/features.h"
#line 1 "/usr/include/features-time64.h"
#line 20 "/usr/include/features-time64.h"
#line 1 "/usr/include/x86_64-linux-gnu/bits/wordsize.h"
#line 21 "/usr/include/features-time64.h"
#line 1 "/usr/include/x86_64-linux-gnu/bits/timesize.h"
#line 19 "/usr/include/x86_64-linux-gnu/bits/timesize.h"
#line 1 "/usr/include/x86_64-linux-gnu/bits/wordsize.h"
#line 20 "/usr/include/x86_64-linux-gnu/bits/timesize.h"
#line 22 "/usr/include/features-time64.h"
#line 393 "/usr/include/features.h"
#line 486 "/usr/include/features.h"
#line 1 "/usr/include/x86_64-linux-gnu/sys/cdefs.h"
#line 559 "/usr/include/x86_64-linux-gnu/sys/cdefs.h"
#line 1 "/usr/include/x86_64-linux-gnu/bits/wordsize.h"
#line 560 "/usr/include/x86_64-linux-gnu/sys/cdefs.h"
#line 1 "/usr/include/x86_64-linux-gnu/bits/long-double.h"
#line 561 "/usr/include/x86_64-linux-gnu/sys/cdefs.h"
#line 487 "/usr/include/features.h"
#line 510 "/usr/include/features.h"
#line 1 "/usr/include/x86_64-linux-gnu/gnu/stubs.h"
#line 10 "/usr/include/x86_64-linux-gnu/gnu/stubs.h"
#line 1 "/usr/include/x86_64-linux-gnu/gnu/stubs-64.h"
#line 11 "/usr/include/x86_64-linux-gnu/gnu/stubs.h"
#line 511 "/usr/include/features.h"
#line 26 "/usr/include/errno.h"


#line 1 "/usr/include/x86_64-linux-gnu/bits/errno.h"
#line 26 "/usr/include/x86_64-linux-gnu/bits/errno.h"
#line 1 "/usr/include/linux/errno.h"
#line 1 "/usr/include/x86_64-linux-gnu/asm/errno.h"
#line 1 "/usr/include/asm-generic/errno.h"




#line 1 "/usr/include/asm-generic/errno-base.h"
#line 6 "/usr/include/asm-generic/errno.h"
#line 2 "/usr/include/x86_64-linux-gnu/asm/errno.h"
#line 2 "/usr/include/linux/errno.h"
#line 27 "/usr/include/x86_64-linux-gnu/bits/errno.h"
#line 29 "/usr/include/errno.h"








extern int *__errno_location (void) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));







extern char *program_invocation_name;
extern char *program_invocation_short_name;

#line 1 "/usr/include/x86_64-linux-gnu/bits/types/error_t.h"
#line 22 "/usr/include/x86_64-linux-gnu/bits/types/error_t.h"
typedef int error_t;
#line 49 "/usr/include/errno.h"




#line 18 "/usr/include/netcdf.h"







typedef int nc_type;
#line 518 "/usr/include/netcdf.h"
 extern const char *
nc_inq_libvers(void);

 extern const char *
nc_strerror(int ncerr);


typedef struct NC_Dispatch NC_Dispatch;
 extern int
nc_def_user_format(int mode_flag, NC_Dispatch *dispatch_table, char *magic_number);

 extern int
nc_inq_user_format(int mode_flag, NC_Dispatch **dispatch_table, char *magic_number);

 extern int
nc__create(const char *path, int cmode, size_t initialsz,
         size_t *chunksizehintp, int *ncidp);

 extern int
nc_create(const char *path, int cmode, int *ncidp);

 extern int
nc__open(const char *path, int mode,
        size_t *chunksizehintp, int *ncidp);

 extern int
nc_open(const char *path, int mode, int *ncidp);


 extern int
nc_inq_path(int ncid, size_t *pathlen, char *path);



 extern int
nc_inq_ncid(int ncid, const char *name, int *grp_ncid);



 extern int
nc_inq_grps(int ncid, int *numgrps, int *ncids);


 extern int
nc_inq_grpname(int ncid, char *name);



 extern int
nc_inq_grpname_full(int ncid, size_t *lenp, char *full_name);


 extern int
nc_inq_grpname_len(int ncid, size_t *lenp);


 extern int
nc_inq_grp_parent(int ncid, int *parent_ncid);


 extern int
nc_inq_grp_ncid(int ncid, const char *grp_name, int *grp_ncid);


 extern int
nc_inq_grp_full_ncid(int ncid, const char *full_name, int *grp_ncid);


 extern int
nc_inq_varids(int ncid, int *nvars, int *varids);



 extern int
nc_inq_dimids(int ncid, int *ndims, int *dimids, int include_parents);



 extern int
nc_inq_typeids(int ncid, int *ntypes, int *typeids);


 extern int
nc_inq_type_equal(int ncid1, nc_type typeid1, int ncid2,
                  nc_type typeid2, int *equal);


 extern int
nc_def_grp(int parent_ncid, const char *name, int *new_ncid);


 extern int
nc_rename_grp(int grpid, const char *name);




 extern int
nc_def_compound(int ncid, size_t size, const char *name, nc_type *typeidp);


 extern int
nc_insert_compound(int ncid, nc_type xtype, const char *name,
                   size_t offset, nc_type field_typeid);


 extern int
nc_insert_array_compound(int ncid, nc_type xtype, const char *name,
                         size_t offset, nc_type field_typeid,
                         int ndims, const int *dim_sizes);


 extern int
nc_inq_type(int ncid, nc_type xtype, char *name, size_t *size);


 extern int
nc_inq_typeid(int ncid, const char *name, nc_type *typeidp);


 extern int
nc_inq_compound(int ncid, nc_type xtype, char *name, size_t *sizep,
                size_t *nfieldsp);


 extern int
nc_inq_compound_name(int ncid, nc_type xtype, char *name);


 extern int
nc_inq_compound_size(int ncid, nc_type xtype, size_t *sizep);


 extern int
nc_inq_compound_nfields(int ncid, nc_type xtype, size_t *nfieldsp);


 extern int
nc_inq_compound_field(int ncid, nc_type xtype, int fieldid, char *name,
                      size_t *offsetp, nc_type *field_typeidp, int *ndimsp,
                      int *dim_sizesp);


 extern int
nc_inq_compound_fieldname(int ncid, nc_type xtype, int fieldid,
                          char *name);


 extern int
nc_inq_compound_fieldindex(int ncid, nc_type xtype, const char *name,
                           int *fieldidp);


 extern int
nc_inq_compound_fieldoffset(int ncid, nc_type xtype, int fieldid,
                            size_t *offsetp);


 extern int
nc_inq_compound_fieldtype(int ncid, nc_type xtype, int fieldid,
                          nc_type *field_typeidp);



 extern int
nc_inq_compound_fieldndims(int ncid, nc_type xtype, int fieldid,
                           int *ndimsp);



 extern int
nc_inq_compound_fielddim_sizes(int ncid, nc_type xtype, int fieldid,
                               int *dim_sizes);


typedef struct {
    size_t len;
    void *p;
} nc_vlen_t;
#line 705 "/usr/include/netcdf.h"
 extern int
nc_def_vlen(int ncid, const char *name, nc_type base_typeid, nc_type *xtypep);


 extern int
nc_inq_vlen(int ncid, nc_type xtype, char *name, size_t *datum_sizep,
            nc_type *base_nc_typep);





 extern int
nc_free_vlen(nc_vlen_t *vl);

 extern int
nc_free_vlens(size_t len, nc_vlen_t vlens[]);


 extern int
nc_put_vlen_element(int ncid, int typeid1, void *vlen_element,
                    size_t len, const void *data);

 extern int
nc_get_vlen_element(int ncid, int typeid1, const void *vlen_element,
                    size_t *len, void *data);





 extern int
nc_free_string(size_t len, char **data);


 extern int
nc_inq_user_type(int ncid, nc_type xtype, char *name, size_t *size,
                 nc_type *base_nc_typep, size_t *nfieldsp, int *classp);


 extern int
nc_put_att(int ncid, int varid, const char *name, nc_type xtype,
           size_t len, const void *op);


 extern int
nc_get_att(int ncid, int varid, const char *name, void *ip);





 extern int
nc_def_enum(int ncid, nc_type base_typeid, const char *name,
            nc_type *typeidp);



 extern int
nc_insert_enum(int ncid, nc_type xtype, const char *name,
               const void *value);



 extern int
nc_inq_enum(int ncid, nc_type xtype, char *name, nc_type *base_nc_typep,
            size_t *base_sizep, size_t *num_membersp);



 extern int
nc_inq_enum_member(int ncid, nc_type xtype, int idx, char *name,
                   void *value);



 extern int
nc_inq_enum_ident(int ncid, nc_type xtype, long long value, char *identifier);




 extern int
nc_def_opaque(int ncid, size_t size, const char *name, nc_type *xtypep);


 extern int
nc_inq_opaque(int ncid, nc_type xtype, char *name, size_t *sizep);


 extern int
nc_put_var(int ncid, int varid, const void *op);


 extern int
nc_get_var(int ncid, int varid, void *ip);


 extern int
nc_put_var1(int ncid, int varid, const size_t *indexp,
            const void *op);


 extern int
nc_get_var1(int ncid, int varid, const size_t *indexp, void *ip);


 extern int
nc_put_vara(int ncid, int varid, const size_t *startp,
            const size_t *countp, const void *op);


 extern int
nc_get_vara(int ncid, int varid, const size_t *startp,
            const size_t *countp, void *ip);


 extern int
nc_put_vars(int ncid, int varid, const size_t *startp,
            const size_t *countp, const ptrdiff_t *stridep,
            const void *op);


 extern int
nc_get_vars(int ncid, int varid, const size_t *startp,
            const size_t *countp, const ptrdiff_t *stridep,
            void *ip);


 extern int
nc_put_varm(int ncid, int varid, const size_t *startp,
            const size_t *countp, const ptrdiff_t *stridep,
            const ptrdiff_t *imapp, const void *op);


 extern int
nc_get_varm(int ncid, int varid, const size_t *startp,
            const size_t *countp, const ptrdiff_t *stridep,
            const ptrdiff_t *imapp, void *ip);





 extern int
nc_def_var_deflate(int ncid, int varid, int shuffle, int deflate,
                   int deflate_level);


 extern int
nc_inq_var_deflate(int ncid, int varid, int *shufflep,
                   int *deflatep, int *deflate_levelp);


 extern int nc_def_var_szip(int ncid, int varid, int options_mask,
                            int pixels_per_block);


 extern int
nc_inq_var_szip(int ncid, int varid, int *options_maskp, int *pixels_per_blockp);



 extern int
nc_def_var_fletcher32(int ncid, int varid, int fletcher32);


 extern int
nc_inq_var_fletcher32(int ncid, int varid, int *fletcher32p);



 extern int
nc_def_var_chunking(int ncid, int varid, int storage, const size_t *chunksizesp);


 extern int
nc_inq_var_chunking(int ncid, int varid, int *storagep, size_t *chunksizesp);



 extern int
nc_def_var_fill(int ncid, int varid, int no_fill, const void *fill_value);


 extern int
nc_inq_var_fill(int ncid, int varid, int *no_fill, void *fill_valuep);


 extern int
nc_def_var_endian(int ncid, int varid, int endian);


 extern int
nc_inq_var_endian(int ncid, int varid, int *endianp);


 extern int
nc_def_var_filter(int ncid, int varid, unsigned int id, size_t nparams, const unsigned int* parms);


 extern int
nc_inq_var_filter(int ncid, int varid, unsigned int* idp, size_t* nparams, unsigned int* params);


 extern int
nc_set_fill(int ncid, int fillmode, int *old_modep);



 extern int
nc_set_default_format(int format, int *old_formatp);


 extern int
nc_set_chunk_cache(size_t size, size_t nelems, float preemption);


 extern int
nc_get_chunk_cache(size_t *sizep, size_t *nelemsp, float *preemptionp);


 extern int
nc_set_var_chunk_cache(int ncid, int varid, size_t size, size_t nelems,
                       float preemption);


 extern int
nc_get_var_chunk_cache(int ncid, int varid, size_t *sizep, size_t *nelemsp,
                       float *preemptionp);

 extern int
nc_redef(int ncid);


 extern int
nc__enddef(int ncid, size_t h_minfree, size_t v_align,
        size_t v_minfree, size_t r_align);

 extern int
nc_enddef(int ncid);

 extern int
nc_sync(int ncid);

 extern int
nc_abort(int ncid);

 extern int
nc_close(int ncid);

 extern int
nc_inq(int ncid, int *ndimsp, int *nvarsp, int *nattsp, int *unlimdimidp);

 extern int
nc_inq_ndims(int ncid, int *ndimsp);

 extern int
nc_inq_nvars(int ncid, int *nvarsp);

 extern int
nc_inq_natts(int ncid, int *nattsp);

 extern int
nc_inq_unlimdim(int ncid, int *unlimdimidp);


 extern int
nc_inq_unlimdims(int ncid, int *nunlimdimsp, int *unlimdimidsp);


 extern int
nc_inq_format(int ncid, int *formatp);


 extern int
nc_inq_format_extended(int ncid, int *formatp, int* modep);



 extern int
nc_def_dim(int ncid, const char *name, size_t len, int *idp);

 extern int
nc_inq_dimid(int ncid, const char *name, int *idp);

 extern int
nc_inq_dim(int ncid, int dimid, char *name, size_t *lenp);

 extern int
nc_inq_dimname(int ncid, int dimid, char *name);

 extern int
nc_inq_dimlen(int ncid, int dimid, size_t *lenp);

 extern int
nc_rename_dim(int ncid, int dimid, const char *name);




 extern int
nc_inq_att(int ncid, int varid, const char *name,
           nc_type *xtypep, size_t *lenp);

 extern int
nc_inq_attid(int ncid, int varid, const char *name, int *idp);

 extern int
nc_inq_atttype(int ncid, int varid, const char *name, nc_type *xtypep);

 extern int
nc_inq_attlen(int ncid, int varid, const char *name, size_t *lenp);

 extern int
nc_inq_attname(int ncid, int varid, int attnum, char *name);

 extern int
nc_copy_att(int ncid_in, int varid_in, const char *name, int ncid_out, int varid_out);

 extern int
nc_rename_att(int ncid, int varid, const char *name, const char *newname);

 extern int
nc_del_att(int ncid, int varid, const char *name);



 extern int
nc_put_att_text(int ncid, int varid, const char *name,
                size_t len, const char *op);

 extern int
nc_get_att_text(int ncid, int varid, const char *name, char *ip);

 extern int
nc_put_att_string(int ncid, int varid, const char *name,
                  size_t len, const char **op);

 extern int
nc_get_att_string(int ncid, int varid, const char *name, char **ip);

 extern int
nc_put_att_uchar(int ncid, int varid, const char *name, nc_type xtype,
                 size_t len, const unsigned char *op);

 extern int
nc_get_att_uchar(int ncid, int varid, const char *name, unsigned char *ip);

 extern int
nc_put_att_schar(int ncid, int varid, const char *name, nc_type xtype,
                 size_t len, const signed char *op);

 extern int
nc_get_att_schar(int ncid, int varid, const char *name, signed char *ip);

 extern int
nc_put_att_short(int ncid, int varid, const char *name, nc_type xtype,
                 size_t len, const short *op);

 extern int
nc_get_att_short(int ncid, int varid, const char *name, short *ip);

 extern int
nc_put_att_int(int ncid, int varid, const char *name, nc_type xtype,
               size_t len, const int *op);

 extern int
nc_get_att_int(int ncid, int varid, const char *name, int *ip);

 extern int
nc_put_att_long(int ncid, int varid, const char *name, nc_type xtype,
                size_t len, const long *op);

 extern int
nc_get_att_long(int ncid, int varid, const char *name, long *ip);

 extern int
nc_put_att_float(int ncid, int varid, const char *name, nc_type xtype,
                 size_t len, const float *op);

 extern int
nc_get_att_float(int ncid, int varid, const char *name, float *ip);

 extern int
nc_put_att_double(int ncid, int varid, const char *name, nc_type xtype,
                  size_t len, const double *op);

 extern int
nc_get_att_double(int ncid, int varid, const char *name, double *ip);

 extern int
nc_put_att_ushort(int ncid, int varid, const char *name, nc_type xtype,
                  size_t len, const unsigned short *op);

 extern int
nc_get_att_ushort(int ncid, int varid, const char *name, unsigned short *ip);

 extern int
nc_put_att_uint(int ncid, int varid, const char *name, nc_type xtype,
                size_t len, const unsigned int *op);

 extern int
nc_get_att_uint(int ncid, int varid, const char *name, unsigned int *ip);

 extern int
nc_put_att_longlong(int ncid, int varid, const char *name, nc_type xtype,
                 size_t len, const long long *op);

 extern int
nc_get_att_longlong(int ncid, int varid, const char *name, long long *ip);

 extern int
nc_put_att_ulonglong(int ncid, int varid, const char *name, nc_type xtype,
                     size_t len, const unsigned long long *op);

 extern int
nc_get_att_ulonglong(int ncid, int varid, const char *name,
                     unsigned long long *ip);





 extern int
nc_def_var(int ncid, const char *name, nc_type xtype, int ndims,
           const int *dimidsp, int *varidp);

 extern int
nc_inq_var(int ncid, int varid, char *name, nc_type *xtypep,
           int *ndimsp, int *dimidsp, int *nattsp);

 extern int
nc_inq_varid(int ncid, const char *name, int *varidp);

 extern int
nc_inq_varname(int ncid, int varid, char *name);

 extern int
nc_inq_vartype(int ncid, int varid, nc_type *xtypep);

 extern int
nc_inq_varndims(int ncid, int varid, int *ndimsp);

 extern int
nc_inq_vardimid(int ncid, int varid, int *dimidsp);

 extern int
nc_inq_varnatts(int ncid, int varid, int *nattsp);

 extern int
nc_rename_var(int ncid, int varid, const char *name);

 extern int
nc_copy_var(int ncid_in, int varid, int ncid_out);
#line 1169 "/usr/include/netcdf.h"
 extern int
nc_put_var1_text(int ncid, int varid, const size_t *indexp, const char *op);

 extern int
nc_get_var1_text(int ncid, int varid, const size_t *indexp, char *ip);

 extern int
nc_put_var1_uchar(int ncid, int varid, const size_t *indexp,
                  const unsigned char *op);

 extern int
nc_get_var1_uchar(int ncid, int varid, const size_t *indexp,
                  unsigned char *ip);

 extern int
nc_put_var1_schar(int ncid, int varid, const size_t *indexp,
                  const signed char *op);

 extern int
nc_get_var1_schar(int ncid, int varid, const size_t *indexp,
                  signed char *ip);

 extern int
nc_put_var1_short(int ncid, int varid, const size_t *indexp,
                  const short *op);

 extern int
nc_get_var1_short(int ncid, int varid, const size_t *indexp,
                  short *ip);

 extern int
nc_put_var1_int(int ncid, int varid, const size_t *indexp, const int *op);

 extern int
nc_get_var1_int(int ncid, int varid, const size_t *indexp, int *ip);

 extern int
nc_put_var1_long(int ncid, int varid, const size_t *indexp, const long *op);

 extern int
nc_get_var1_long(int ncid, int varid, const size_t *indexp, long *ip);

 extern int
nc_put_var1_float(int ncid, int varid, const size_t *indexp, const float *op);

 extern int
nc_get_var1_float(int ncid, int varid, const size_t *indexp, float *ip);

 extern int
nc_put_var1_double(int ncid, int varid, const size_t *indexp, const double *op);

 extern int
nc_get_var1_double(int ncid, int varid, const size_t *indexp, double *ip);

 extern int
nc_put_var1_ushort(int ncid, int varid, const size_t *indexp,
                   const unsigned short *op);

 extern int
nc_get_var1_ushort(int ncid, int varid, const size_t *indexp,
                   unsigned short *ip);

 extern int
nc_put_var1_uint(int ncid, int varid, const size_t *indexp,
                 const unsigned int *op);

 extern int
nc_get_var1_uint(int ncid, int varid, const size_t *indexp,
                 unsigned int *ip);

 extern int
nc_put_var1_longlong(int ncid, int varid, const size_t *indexp,
                     const long long *op);

 extern int
nc_get_var1_longlong(int ncid, int varid, const size_t *indexp,
                  long long *ip);

 extern int
nc_put_var1_ulonglong(int ncid, int varid, const size_t *indexp,
                   const unsigned long long *op);

 extern int
nc_get_var1_ulonglong(int ncid, int varid, const size_t *indexp,
                   unsigned long long *ip);

 extern int
nc_put_var1_string(int ncid, int varid, const size_t *indexp,
                   const char **op);

 extern int
nc_get_var1_string(int ncid, int varid, const size_t *indexp,
                   char **ip);




 extern int
nc_put_vara_text(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const char *op);

 extern int
nc_get_vara_text(int ncid, int varid, const size_t *startp,
                 const size_t *countp, char *ip);

 extern int
nc_put_vara_uchar(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const unsigned char *op);

 extern int
nc_get_vara_uchar(int ncid, int varid, const size_t *startp,
                  const size_t *countp, unsigned char *ip);

 extern int
nc_put_vara_schar(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const signed char *op);

 extern int
nc_get_vara_schar(int ncid, int varid, const size_t *startp,
                  const size_t *countp, signed char *ip);

 extern int
nc_put_vara_short(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const short *op);

 extern int
nc_get_vara_short(int ncid, int varid, const size_t *startp,
                  const size_t *countp, short *ip);

 extern int
nc_put_vara_int(int ncid, int varid, const size_t *startp,
                const size_t *countp, const int *op);

 extern int
nc_get_vara_int(int ncid, int varid, const size_t *startp,
                const size_t *countp, int *ip);

 extern int
nc_put_vara_long(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const long *op);

 extern int
nc_get_vara_long(int ncid, int varid,
        const size_t *startp, const size_t *countp, long *ip);

 extern int
nc_put_vara_float(int ncid, int varid,
        const size_t *startp, const size_t *countp, const float *op);

 extern int
nc_get_vara_float(int ncid, int varid,
        const size_t *startp, const size_t *countp, float *ip);

 extern int
nc_put_vara_double(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const double *op);

 extern int
nc_get_vara_double(int ncid, int varid, const size_t *startp,
                   const size_t *countp, double *ip);

 extern int
nc_put_vara_ushort(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const unsigned short *op);

 extern int
nc_get_vara_ushort(int ncid, int varid, const size_t *startp,
                   const size_t *countp, unsigned short *ip);

 extern int
nc_put_vara_uint(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const unsigned int *op);

 extern int
nc_get_vara_uint(int ncid, int varid, const size_t *startp,
                 const size_t *countp, unsigned int *ip);

 extern int
nc_put_vara_longlong(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const long long *op);

 extern int
nc_get_vara_longlong(int ncid, int varid, const size_t *startp,
                  const size_t *countp, long long *ip);

 extern int
nc_put_vara_ulonglong(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const unsigned long long *op);

 extern int
nc_get_vara_ulonglong(int ncid, int varid, const size_t *startp,
                   const size_t *countp, unsigned long long *ip);

 extern int
nc_put_vara_string(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const char **op);

 extern int
nc_get_vara_string(int ncid, int varid, const size_t *startp,
                   const size_t *countp, char **ip);




 extern int
nc_put_vars_text(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        const char *op);

 extern int
nc_get_vars_text(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        char *ip);

 extern int
nc_put_vars_uchar(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        const unsigned char *op);

 extern int
nc_get_vars_uchar(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        unsigned char *ip);

 extern int
nc_put_vars_schar(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        const signed char *op);

 extern int
nc_get_vars_schar(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        signed char *ip);

 extern int
nc_put_vars_short(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        const short *op);

 extern int
nc_get_vars_short(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  short *ip);

 extern int
nc_put_vars_int(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        const int *op);

 extern int
nc_get_vars_int(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        int *ip);

 extern int
nc_put_vars_long(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        const long *op);

 extern int
nc_get_vars_long(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        long *ip);

 extern int
nc_put_vars_float(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        const float *op);

 extern int
nc_get_vars_float(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        float *ip);

 extern int
nc_put_vars_double(int ncid, int varid,
        const size_t *startp, const size_t *countp, const ptrdiff_t *stridep,
        const double *op);

 extern int
nc_get_vars_double(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   double *ip);

 extern int
nc_put_vars_ushort(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const unsigned short *op);

 extern int
nc_get_vars_ushort(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   unsigned short *ip);

 extern int
nc_put_vars_uint(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const ptrdiff_t *stridep,
                 const unsigned int *op);

 extern int
nc_get_vars_uint(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const ptrdiff_t *stridep,
                 unsigned int *ip);

 extern int
nc_put_vars_longlong(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const long long *op);

 extern int
nc_get_vars_longlong(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  long long *ip);

 extern int
nc_put_vars_ulonglong(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const unsigned long long *op);

 extern int
nc_get_vars_ulonglong(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   unsigned long long *ip);

 extern int
nc_put_vars_string(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const char **op);

 extern int
nc_get_vars_string(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   char **ip);




 extern int
nc_put_varm_text(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const ptrdiff_t *stridep,
                 const ptrdiff_t *imapp, const char *op);

 extern int
nc_get_varm_text(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const ptrdiff_t *stridep,
                 const ptrdiff_t *imapp, char *ip);

 extern int
nc_put_varm_uchar(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t *imapp, const unsigned char *op);

 extern int
nc_get_varm_uchar(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t *imapp, unsigned char *ip);

 extern int
nc_put_varm_schar(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t *imapp, const signed char *op);

 extern int
nc_get_varm_schar(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t *imapp, signed char *ip);

 extern int
nc_put_varm_short(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t *imapp, const short *op);

 extern int
nc_get_varm_short(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t *imapp, short *ip);

 extern int
nc_put_varm_int(int ncid, int varid, const size_t *startp,
                const size_t *countp, const ptrdiff_t *stridep,
                const ptrdiff_t *imapp, const int *op);

 extern int
nc_get_varm_int(int ncid, int varid, const size_t *startp,
                const size_t *countp, const ptrdiff_t *stridep,
                const ptrdiff_t *imapp, int *ip);

 extern int
nc_put_varm_long(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const ptrdiff_t *stridep,
                 const ptrdiff_t *imapp, const long *op);

 extern int
nc_get_varm_long(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const ptrdiff_t *stridep,
                 const ptrdiff_t *imapp, long *ip);

 extern int
nc_put_varm_float(int ncid, int varid,const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t *imapp, const float *op);

 extern int
nc_get_varm_float(int ncid, int varid,const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t *imapp, float *ip);

 extern int
nc_put_varm_double(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const ptrdiff_t *imapp, const double *op);

 extern int
nc_get_varm_double(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const ptrdiff_t * imapp, double *ip);

 extern int
nc_put_varm_ushort(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const ptrdiff_t * imapp, const unsigned short *op);

 extern int
nc_get_varm_ushort(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const ptrdiff_t * imapp, unsigned short *ip);

 extern int
nc_put_varm_uint(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const ptrdiff_t *stridep,
                 const ptrdiff_t * imapp, const unsigned int *op);

 extern int
nc_get_varm_uint(int ncid, int varid, const size_t *startp,
                 const size_t *countp, const ptrdiff_t *stridep,
                 const ptrdiff_t * imapp, unsigned int *ip);

 extern int
nc_put_varm_longlong(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t * imapp, const long long *op);

 extern int
nc_get_varm_longlong(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t * imapp, long long *ip);

 extern int
nc_put_varm_ulonglong(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const ptrdiff_t * imapp, const unsigned long long *op);

 extern int
nc_get_varm_ulonglong(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const ptrdiff_t * imapp, unsigned long long *ip);

 extern int
nc_put_varm_string(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const ptrdiff_t * imapp, const char **op);

 extern int
nc_get_varm_string(int ncid, int varid, const size_t *startp,
                   const size_t *countp, const ptrdiff_t *stridep,
                   const ptrdiff_t * imapp, char **ip);




 extern int
nc_put_var_text(int ncid, int varid, const char *op);

 extern int
nc_get_var_text(int ncid, int varid, char *ip);

 extern int
nc_put_var_uchar(int ncid, int varid, const unsigned char *op);

 extern int
nc_get_var_uchar(int ncid, int varid, unsigned char *ip);

 extern int
nc_put_var_schar(int ncid, int varid, const signed char *op);

 extern int
nc_get_var_schar(int ncid, int varid, signed char *ip);

 extern int
nc_put_var_short(int ncid, int varid, const short *op);

 extern int
nc_get_var_short(int ncid, int varid, short *ip);

 extern int
nc_put_var_int(int ncid, int varid, const int *op);

 extern int
nc_get_var_int(int ncid, int varid, int *ip);

 extern int
nc_put_var_long(int ncid, int varid, const long *op);

 extern int
nc_get_var_long(int ncid, int varid, long *ip);

 extern int
nc_put_var_float(int ncid, int varid, const float *op);

 extern int
nc_get_var_float(int ncid, int varid, float *ip);

 extern int
nc_put_var_double(int ncid, int varid, const double *op);

 extern int
nc_get_var_double(int ncid, int varid, double *ip);

 extern int
nc_put_var_ushort(int ncid, int varid, const unsigned short *op);

 extern int
nc_get_var_ushort(int ncid, int varid, unsigned short *ip);

 extern int
nc_put_var_uint(int ncid, int varid, const unsigned int *op);

 extern int
nc_get_var_uint(int ncid, int varid, unsigned int *ip);

 extern int
nc_put_var_longlong(int ncid, int varid, const long long *op);

 extern int
nc_get_var_longlong(int ncid, int varid, long long *ip);

 extern int
nc_put_var_ulonglong(int ncid, int varid, const unsigned long long *op);

 extern int
nc_get_var_ulonglong(int ncid, int varid, unsigned long long *ip);

 extern int
nc_put_var_string(int ncid, int varid, const char **op);

 extern int
nc_get_var_string(int ncid, int varid, char **ip);


 extern int
nc_put_att_ubyte(int ncid, int varid, const char *name, nc_type xtype,
                 size_t len, const unsigned char *op);
 extern int
nc_get_att_ubyte(int ncid, int varid, const char *name,
                 unsigned char *ip);
 extern int
nc_put_var1_ubyte(int ncid, int varid, const size_t *indexp,
                  const unsigned char *op);
 extern int
nc_get_var1_ubyte(int ncid, int varid, const size_t *indexp,
                  unsigned char *ip);
 extern int
nc_put_vara_ubyte(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const unsigned char *op);
 extern int
nc_get_vara_ubyte(int ncid, int varid, const size_t *startp,
                  const size_t *countp, unsigned char *ip);
 extern int
nc_put_vars_ubyte(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const unsigned char *op);
 extern int
nc_get_vars_ubyte(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  unsigned char *ip);
 extern int
nc_put_varm_ubyte(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t * imapp, const unsigned char *op);
 extern int
nc_get_varm_ubyte(int ncid, int varid, const size_t *startp,
                  const size_t *countp, const ptrdiff_t *stridep,
                  const ptrdiff_t * imapp, unsigned char *ip);
 extern int
nc_put_var_ubyte(int ncid, int varid, const unsigned char *op);
 extern int
nc_get_var_ubyte(int ncid, int varid, unsigned char *ip);




 extern int
nc_set_log_level(int new_level);






 extern int
nc_show_metadata(int ncid);




 extern int
nc_delete(const char *path);
#line 1786 "/usr/include/netcdf.h"
 extern int
nc__create_mp(const char *path, int cmode, size_t initialsz, int basepe,
         size_t *chunksizehintp, int *ncidp);

 extern int
nc__open_mp(const char *path, int mode, int basepe,
        size_t *chunksizehintp, int *ncidp);

 extern int
nc_delete_mp(const char *path, int basepe);

 extern int
nc_set_base_pe(int ncid, int pe);

 extern int
nc_inq_base_pe(int ncid, int *pe);


 extern int
nctypelen(nc_type datatype);
#line 1829 "/usr/include/netcdf.h"
 extern int ncerr;
#line 1843 "/usr/include/netcdf.h"
 extern int ncopts;

 extern void
nc_advise(const char *cdf_routine_name, int err, const char *fmt,...);






typedef int nclong;

 extern int
nccreate(const char* path, int cmode);

 extern int
ncopen(const char* path, int mode);

 extern int
ncsetfill(int ncid, int fillmode);

 extern int
ncredef(int ncid);

 extern int
ncendef(int ncid);

 extern int
ncsync(int ncid);

 extern int
ncabort(int ncid);

 extern int
ncclose(int ncid);

 extern int
ncinquire(int ncid, int *ndimsp, int *nvarsp, int *nattsp, int *unlimdimp);

 extern int
ncdimdef(int ncid, const char *name, long len);

 extern int
ncdimid(int ncid, const char *name);

 extern int
ncdiminq(int ncid, int dimid, char *name, long *lenp);

 extern int
ncdimrename(int ncid, int dimid, const char *name);

 extern int
ncattput(int ncid, int varid, const char *name, nc_type xtype,
        int len, const void *op);

 extern int
ncattinq(int ncid, int varid, const char *name, nc_type *xtypep, int *lenp);

 extern int
ncattget(int ncid, int varid, const char *name, void *ip);

 extern int
ncattcopy(int ncid_in, int varid_in, const char *name, int ncid_out,
        int varid_out);

 extern int
ncattname(int ncid, int varid, int attnum, char *name);

 extern int
ncattrename(int ncid, int varid, const char *name, const char *newname);

 extern int
ncattdel(int ncid, int varid, const char *name);

 extern int
ncvardef(int ncid, const char *name, nc_type xtype,
        int ndims, const int *dimidsp);

 extern int
ncvarid(int ncid, const char *name);

 extern int
ncvarinq(int ncid, int varid, char *name, nc_type *xtypep,
        int *ndimsp, int *dimidsp, int *nattsp);

 extern int
ncvarput1(int ncid, int varid, const long *indexp, const void *op);

 extern int
ncvarget1(int ncid, int varid, const long *indexp, void *ip);

 extern int
ncvarput(int ncid, int varid, const long *startp, const long *countp,
        const void *op);

 extern int
ncvarget(int ncid, int varid, const long *startp, const long *countp,
        void *ip);

 extern int
ncvarputs(int ncid, int varid, const long *startp, const long *countp,
        const long *stridep, const void *op);

 extern int
ncvargets(int ncid, int varid, const long *startp, const long *countp,
        const long *stridep, void *ip);

 extern int
ncvarputg(int ncid, int varid, const long *startp, const long *countp,
        const long *stridep, const long *imapp, const void *op);

 extern int
ncvargetg(int ncid, int varid, const long *startp, const long *countp,
        const long *stridep, const long *imapp, void *ip);

 extern int
ncvarrename(int ncid, int varid, const char *name);

 extern int
ncrecinq(int ncid, int *nrecvarsp, int *recvaridsp, long *recsizesp);

 extern int
ncrecget(int ncid, long recnum, void **datap);

 extern int
ncrecput(int ncid, long recnum, void *const *datap);




 extern int nc_initialize(void);





 extern int nc_finalize(void);
#line 8 "./netcdf_vertex_bas.h"

#line 25 "./netcdf_vertex_bas.h"

#line 25 "./netcdf_vertex_bas.h"
int nc_err;



scalar * scalar_list_nc;
char file_nc[80];


int ncid;
int t_varid;


int nc_varid[1000];
char * nc_varname[1000];
int nvarout = 0;
int nc_rec = -1;


int nl_tmp = 1;

void create_nc()
{


   int x_dimid, y_dimid, lvl_dimid, rec_dimid;
   int y_varid, x_varid;
   int dimids[3];


   if ((nc_err = nc_create(file_nc, 
#line 54 "./netcdf_vertex_bas.h"
                                   0x0000
#line 54 "./netcdf_vertex_bas.h"
                                             , &ncid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};






   if ((nc_err = nc_def_dim(ncid, "y", N+1, &y_dimid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
   if ((nc_err = nc_def_dim(ncid, "x", N+1, &x_dimid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
   if ((nc_err = nc_def_dim(ncid, "time", 
#line 66 "./netcdf_vertex_bas.h"
                                           0L
#line 66 "./netcdf_vertex_bas.h"
                                                       , &rec_dimid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};







   if ((nc_err = nc_def_var(ncid, "time", 
#line 75 "./netcdf_vertex_bas.h"
                                           5
#line 75 "./netcdf_vertex_bas.h"
                                                   , 1, &rec_dimid,
              &t_varid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
   if ((nc_err = nc_def_var(ncid, "y", 
#line 78 "./netcdf_vertex_bas.h"
                                         5
#line 78 "./netcdf_vertex_bas.h"
                                                 , 1, &y_dimid,
              &y_varid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
   if ((nc_err = nc_def_var(ncid, "x", 
#line 81 "./netcdf_vertex_bas.h"
                                         5
#line 81 "./netcdf_vertex_bas.h"
                                                 , 1, &x_dimid,
              &x_varid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};





   dimids[0] = rec_dimid;



   dimids[1] = y_dimid;
   dimids[2] = x_dimid;



   {scalar*_i=(scalar*)( scalar_list_nc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){{

       if ((nc_err = nc_def_var(ncid, _attribute[s.i].name, 
#line 100 "./netcdf_vertex_bas.h"
                                             5
#line 100 "./netcdf_vertex_bas.h"
                                                     , 3,
                                dimids, &nc_varid[nvarout])))
         {printf("Error: %s\n", nc_strerror(nc_err)); return;};

       nvarout += 1;


   }}}
#line 118 "./netcdf_vertex_bas.h"
   if ((nc_err = nc_enddef(ncid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};


   float yc[N+1], xc[N+1];
   double Delta = L0*1.0/N;
   for (int i = 0; i < N+1; i++){
      yc[i] = Y0 + i*Delta;
      xc[i] = X0 + i*Delta;
   }

   if ((nc_err = nc_put_var_float(ncid, y_varid, &yc[0])))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
   if ((nc_err = nc_put_var_float(ncid, x_varid, &xc[0])))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};


   if ((nc_err = nc_close(ncid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
   printf("*** SUCCESS creating example file %s!\n", file_nc);

}


struct OutputNetcdf {
  int it;
  int n;
  bool linear;
} OutputNetcdf;

void write_nc(struct OutputNetcdf p) {
  if (p.n == 0) p.n = N + 1;

  if (pid() == 0) {

    if ((nc_err = nc_open(file_nc, 
#line 153 "./netcdf_vertex_bas.h"
                                  0x0001
#line 153 "./netcdf_vertex_bas.h"
                                          , &ncid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
  }


  nc_rec += 1;
  float loctime = t;

  size_t startt[1], countt[1];
  startt[0] = nc_rec;
  countt[0] = 1;
  if (pid() == 0) {
    if ((nc_err = nc_put_vara_float(ncid, t_varid, startt, countt,
                                    &loctime)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
  }



  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));



  size_t start[3], count[3];
#line 192 "./netcdf_vertex_bas.h"
  start[0] = nc_rec;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = p.n;
  count[2] = p.n;


  int nv = -1;


  {scalar*_i=(scalar*)( scalar_list_nc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){{
    nv += 1;

    for (int j = 0; j < p.n; j++) {
      for (int i = 0; i < p.n; i++) {
        field[j][i] = HUGE;
      }
    }

    {foreach_vertex(){

      field[(point.j - 2)][(point.i - 2)] = val(s,0,0,0);
    }end_foreach_vertex();}
#line 236 "./netcdf_vertex_bas.h"
    if (pid() == 0) {
#if _MPI
        MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif
#line 248 "./netcdf_vertex_bas.h"
     if ((nc_err = nc_put_vara_float(ncid, nc_varid[nv], start, count,
                 &field[0][0])))
         {printf("Error: %s\n", nc_strerror(nc_err)); return;};



  }
#if _MPI
  else
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif

  }}}
  matrix_free (field);



  if (pid() == 0) {
    if ((nc_err = nc_close(ncid)))
      {printf("Error: %s\n", nc_strerror(nc_err)); return;};
  }

}






void read_nc(scalar * list_in, char* file_in){

  int i, ret;
  int ncfile, ndims, nvars, ngatts, unlimited;
  int var_ndims, var_natts;
  nc_type type;
  char varname[
#line 283 "./netcdf_vertex_bas.h"
              256
#line 283 "./netcdf_vertex_bas.h"
                         +1];
  int *dimids=NULL;

  int Nloc = N+1;
  float ** field = matrix_new (Nloc, Nloc, sizeof(float));

  if ((nc_err = nc_open(file_in, 
#line 289 "./netcdf_vertex_bas.h"
                                0x0000
#line 289 "./netcdf_vertex_bas.h"
                                          , &ncfile)))
    {printf("Error: %s\n", nc_strerror(nc_err)); return;};

  if ((nc_err = nc_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited)))
    {printf("Error: %s\n", nc_strerror(nc_err)); return;};

  size_t start[3], count[3];
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = Nloc;
  count[2] = Nloc;


  {scalar*_i=(scalar*)( list_in);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){{
    for(i=0; i<nvars; i++) {


  if ((nc_err = nc_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                          &var_natts)))
    {printf("Error: %s\n", nc_strerror(nc_err)); return;};


      if (strcmp(varname,_attribute[s.i].name) == 0) {
        fprintf(fout,"Reading variable  %s!\n", _attribute[s.i].name);

          if ((nc_err = nc_get_vara_float(ncfile, i, start, count,
                                                 &field[0][0])))
            {printf("Error: %s\n", nc_strerror(nc_err)); return;};

          {foreach_vertex(){
            val(s,0,0,0) = field[(point.j - 2)][(point.i - 2)];
          }end_foreach_vertex();}

        }


    }
  }}}

  matrix_free (field);

  if ((nc_err = nc_close(ncfile)))
    {printf("Error: %s\n", nc_strerror(nc_err)); return;};

  boundary_internal ((scalar *)list_in, "./netcdf_vertex_bas.h", 336);

}

static int cleanup_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 1234567890);*ip=i;*tp=t;return ret;}      static int cleanup_0(const int i,const double t,Event *_ev){tracing("cleanup_0","./netcdf_vertex_bas.h",340);
{
  pfree(scalar_list_nc,__func__,__FILE__,__LINE__);
}{end_tracing("cleanup_0","./netcdf_vertex_bas.h",343);return 0;}end_tracing("cleanup_0","./netcdf_vertex_bas.h",343);}
#line 24 "qg.c"

char* fileout = "vars.nc";

int main(int argc,char* argv[]) {_init_solver();


  if (argc == 2) {
    read_params(argv[1]);
  } else {
    read_params("params.in");
  }

  if (sbc == -1) {
    periodic(right);
    periodic(top);
  }

  create_outdir();

  init_grid (N);
  size(L0);

  run();
free_solver();}





static int init_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int init_0(const int i,const double t,Event *_ev){tracing("init_0","qg.c",53); {

  foreach_vertex_stencil()
    {_stencil_val_a(psi,0,0,0);  }end_foreach_vertex_stencil();

  {
#line 55
foreach_vertex()
    val(psi,0,0,0) = 0;end_foreach_vertex();}
  FILE * fp;
  if ((fp = fopen("restart.nc", "r"))) {
    read_nc(((scalar[]){psi,{-1}}), "restart.nc");
    fclose(fp);
  }

  boundary_internal ((scalar *)((scalar[]){psi,{-1}}), "qg.c", 63);

}{end_tracing("init_0","qg.c",65);return 0;}end_tracing("init_0","qg.c",65);}




static int write_const_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 0);*ip=i;*tp=t;return ret;}      static int write_const(const int i,const double t,Event *_ev){tracing("write_const","qg.c",70); {
  backup_config();

  sprintf (file_nc,"%s%s", dpath, fileout);
  scalar_list_nc = list_copy(((scalar[]){psi, q,{-1}}));
  create_nc();
}{end_tracing("write_const","qg.c",76);return 0;}end_tracing("write_const","qg.c",76);}

static int writestdout_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int writestdout(const int i,const double t,Event *_ev){tracing("writestdout","qg.c",78); {
  double ke = 0;
  foreach_vertex_stencil()
    {_stencil_val(psi,0,0,0);_stencil_val(psi,1,0,0); _stencil_val(psi,-1,0,0); _stencil_val(psi,0,1,0); _stencil_val(psi,0,-1,0);_stencil_val(psi,0,0,0);       }end_foreach_vertex_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:ke)){
#line 80
foreach_vertex()
    ke -= 0.5*val(psi,0,0,0)*(val(psi,1,0,0) + val(psi,-1,0,0) + val(psi,0,1,0) + val(psi,0,-1,0) - 4*val(psi,0,0,0))/(sq(Delta))*sq(Delta);end_foreach_vertex();mpi_all_reduce_array(&ke,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

  
#line 83
fprintf (fout,"i = %i, dt = %g, t = %g, ke_1 = %g\n", i, dt, t, ke);
}{end_tracing("writestdout","qg.c",84);return 0;}end_tracing("writestdout","qg.c",84);}

static int output_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=( i += 100);*ip=i;*tp=t;return ret;}static int output_expr1(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=( i <= iend);*ip=i;*tp=t;return ret;}static int output_expr2(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int output(const int i,const double t,Event *_ev){tracing("output","qg.c",86);{
  fprintf(fout,"write file\n");
  write_nc((struct OutputNetcdf){0});
}{end_tracing("output","qg.c",89);return 0;}end_tracing("output","qg.c",89);}
#line 2 "ast/init_solver.h"

static void _init_solver (void)
{
  void init_solver();
  datasize=3*sizeof(double);init_solver();
  multigrid_methods();{  event_register((Event){0,1,defaults,{defaults_expr0},((int *)0),((double *)0),"/home/lennard/basilisk/src/predictor-corrector.h",34,"defaults"});
  event_register((Event){0,1,defaults_0,{defaults_0_expr0},((int *)0),((double *)0),"./qg.h",323,"defaults"});
  event_register((Event){0,1,init,{init_expr0},((int *)0),((double *)0),"./qg.h",337,"init"});
  event_register((Event){0,1,init_0,{init_0_expr0},((int *)0),((double *)0),"qg.c",53,"init"});
  event_register((Event){0,1,write_const,{write_const_expr0},((int *)0),((double *)0),"qg.c",70,"write_const"});
  event_register((Event){0,1,writestdout,{writestdout_expr0},((int *)0),((double *)0),"qg.c",78,"writestdout"});
  event_register((Event){0,3,output,{output_expr0,output_expr1,output_expr2},((int *)0),((double *)0),"qg.c",86,"output"});


    

    init_const_vector((vector){{_NVARMAX+0},{_NVARMAX+1}},"zerof",(double[]) {0.,0.,0.});
  init_const_vector((vector){{_NVARMAX+2},{_NVARMAX+3}},"unityf",(double[]) {1.,1.,1.});
  init_const_scalar((scalar){_NVARMAX+4},"unity", 1.);
  init_const_scalar((scalar){_NVARMAX+5},"zeroc", 0.);
  init_vertex_scalar((scalar){0},"psi");
  init_vertex_scalar((scalar){1},"q");
  init_vertex_scalar((scalar){2},"zeta");
  event_register((Event){0,1,cleanup,{cleanup_expr0},((int *)0),((double *)0),"./qg.h",351,"cleanup"});
  event_register((Event){0,1,cleanup_0,{cleanup_0_expr0},((int *)0),((double *)0),"./netcdf_vertex_bas.h",340,"cleanup"});

#line 13
}
  set_fpe();
}
