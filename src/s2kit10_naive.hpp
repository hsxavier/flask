#ifndef S2KIT10NAIVE_H
#define S2KIT10NAIVE_H 1

void ArcCosEvalPts(int n, double *eval_pts);
void EvalPts(int n, double *eval_pts);
void vec_add(double *data1, double *data2, double *result, int n);
void vec_mul(double scalar, double *data1, double *result, int n);
void vec_pt_mul(double *data1, double *data2, double *result, int n);
double L2_an(int m, int l);
double L2_cn(int m, int l);
void Pmm_L2(int m, double *eval_pts, int n, double *result);
void PmlTableGen(int bw, int m, double *storeplm, double *workspace);
void Naive_SynthesizeX(double *coeffs, int bw, int m, double *result, double *plmtable);

void makeweights(int bw, double *weights);
void Naive_AnalysisX(double *data, int bw, int m, double *weights, double *result, double *plmtable, double *workspace);

double suppress(double l, double lsup, double supindex);
double *GetCl4DLT(double *Clin, double *ll, int Clinsize, double lsup, double supindex, int lmax);
void GetAllLs(double *ll, double *Clin, int Clinsize, double *Clout, int lmax, int extrapol=0);
void ModCl4DLT(double *Clout, int lmax, double lsup, double supindex);
double unsuppress(double l, double lsup, double supindex);
void ApplyClFactors(double *Cl, int ClLength, double lsup=-1, double supindex=0);
 
#endif
