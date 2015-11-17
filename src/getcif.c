#include <R.h>
#include "methas.h"

void fexitc(const char *msg);

extern Cifns AreaIntCifns, BadGeyCifns, DgsCifns, DiggraCifns, 
  FikselCifns, GeyerCifns, HardcoreCifns, 
  LennardCifns, LookupCifns, 
  SoftcoreCifns, StraussCifns, StraussHardCifns, 
  MultiStraussCifns, MultiStraussHardCifns, MultiHardCifns,
  TripletsCifns, PenttinenCifns;

Cifns NullCifns = NULL_CIFNS;

typedef struct CifPair {
  char *name;
  Cifns *p;
} CifPair;

CifPair CifTable[] = { 
  {"areaint",   &AreaIntCifns},
  {"badgey",    &BadGeyCifns},
  {"dgs",       &DgsCifns},
  {"diggra",    &DiggraCifns},
  {"geyer",     &GeyerCifns},
  {"fiksel",    &FikselCifns},
  {"hardcore",  &HardcoreCifns},
  {"lookup",    &LookupCifns},
  {"lennard",   &LennardCifns},
  {"multihard", &MultiHardCifns},
  {"penttinen", &PenttinenCifns},
  {"sftcr",     &SoftcoreCifns},
  {"strauss",   &StraussCifns},
  {"straush",   &StraussHardCifns},
  {"straussm",  &MultiStraussCifns},
  {"straushm",  &MultiStraussHardCifns},
  {"triplets",  &TripletsCifns},
  {(char *) NULL, (Cifns *) NULL}
};

Cifns getcif(cifname) 
     char *cifname;
{
  int i;
  CifPair cp;
  for(i = 0; CifTable[i].name; i++) {
    cp = CifTable[i];
    if(strcmp(cifname, cp.name) == 0)
      return(*(cp.p));
  }
  fexitc("Unrecognised cif name; bailing out.\n");
  /* control never passes to here, but compilers don't know that */
  return(NullCifns);
}

/* R interface function, to check directly whether cif is recognised */

void knownCif(cifname, answer) 
     char** cifname;
     int* answer;
{
  int i;
  CifPair cp;
  for(i = 0; CifTable[i].name; i++) {
    cp = CifTable[i];
    if(strcmp(*cifname, cp.name) == 0) {
      *answer = 1;
      return;
    }
  }
  *answer = 0;
  return;
}

