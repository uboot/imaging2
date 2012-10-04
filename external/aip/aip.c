// http://www.cap-lore.com/MathPhys/IP/
// see also: http://www.cap-lore.com/MathPhys/IP/


#include <stdio.h>
typedef double real; // could be double.
#define C const
#define bigReal 1.E38  // FLT_MAX, DBL_MAX if double above
#include <limits.h>
typedef struct{real x; real y;} point;
real inter(point * a, int na, point * b, int nb){
  typedef struct{point min; point max;} box;
  typedef long long hp;
  typedef struct{long x; long y;} ipoint;
  typedef struct{long mn; long mx;} rng;
  typedef struct{ipoint ip; rng rx; rng ry; short in;} vertex;
  vertex ipa[na+1], ipb[nb+1];
  box B = {{bigReal, bigReal}, {-bigReal, -bigReal}};
  double ascale;
  void range(point * x, int c){
    while(c--){
      void bd(real * X, real y){*X = *X<y ? *X:y;}
      void bu(real * X, real y){*X = *X>y ? *X:y;}
      
      bd(&B.min.x, x[c].x); bu(&B.max.x, x[c].x);
      bd(&B.min.y, x[c].y); bu(&B.max.y, x[c].y);}}
  if(na < 3 || nb < 3) return 0;
  range(a, na); range(b, nb);
  {  const real gamut = 500000000., mid = gamut/2.;
     real rngx = B.max.x - B.min.x, sclx = gamut/rngx,
          rngy = B.max.y - B.min.y, scly = gamut/rngy;
     {void fit(point * x, int cx, vertex * ix, int fudge){
        {int c=cx; while(c--){
          ix[c].ip.x = (long)((x[c].x - B.min.x)*sclx - mid)&~7|fudge|c&1;
          ix[c].ip.y = (long)((x[c].y - B.min.y)*scly - mid)&~7|fudge;
          }}
        ix[0].ip.y += cx&1;
        ix[cx] = ix[0];
        {int c=cx; while(c--) {
           ix[c].rx = ix[c].ip.x < ix[c+1].ip.x ?
             (rng){ix[c].ip.x,ix[c+1].ip.x}:(rng){ix[c+1].ip.x,ix[c].ip.x};
           ix[c].ry = ix[c].ip.y < ix[c+1].ip.y ?
             (rng){ix[c].ip.y,ix[c+1].ip.y}:(rng){ix[c+1].ip.y,ix[c].ip.y};
           ix[c].in=0;}}}
      fit(a, na, ipa, 0); fit(b, nb, ipb, 2);}
  ascale = sclx*scly;}
  { hp area(ipoint a, ipoint p, ipoint q)
        {return (hp)p.x*q.y - (hp)p.y*q.x +
            (hp)a.x*(p.y - q.y) + (hp)a.y*(q.x - p.x);}
    hp s = 0; int j, k;
    void cntrib(ipoint f, ipoint t, short w)
       {s += (hp)w*(t.x-f.x)*(t.y+f.y)/2;}
    int ovl(rng p, rng q){return p.mn < q.mx && q.mn < p.mx;}
    for(j=0; j<na; ++j) for(k=0; k<nb; ++k)
      if(ovl(ipa[j].rx, ipb[k].rx) && ovl(ipa[j].ry, ipb[k].ry)){
         hp a1 = -area(ipa[j].ip, ipb[k].ip, ipb[k+1].ip),
           a2 = area(ipa[j+1].ip, ipb[k].ip, ipb[k+1].ip);
        {int o = a1<0; if(o == a2<0)
          {hp a3 = area(ipb[k].ip, ipa[j].ip, ipa[j+1].ip),
            a4 = -area(ipb[k+1].ip, ipa[j].ip, ipa[j+1].ip);
          if(a3<0 == a4<0) {
            void cross(vertex * a, vertex * b, vertex * c, vertex * d,
               double a1, double a2, double a3, double a4){
               real r1=a1/((real)a1+a2), r2 = a3/((real)a3+a4);
               cntrib((ipoint){a->ip.x + r1*(b->ip.x-a->ip.x), 
                            a->ip.y + r1*(b->ip.y-a->ip.y)},
                                b->ip, 1);
               cntrib(d->ip, (ipoint){
                 c->ip.x + r2*(d->ip.x - c->ip.x), 
                 c->ip.y + r2*(d->ip.y - c->ip.y)}, 1);
                ++a->in; --c->in;}
            if(o) cross(&ipa[j], &ipa[j+1], &ipb[k], &ipb[k+1], a1, a2, a3, a4);
            else cross(&ipb[k], &ipb[k+1], &ipa[j], &ipa[j+1], a3, a4, a1, a2);
   }}}}
   {void inness(vertex * P, int cP, vertex * Q, int cQ){
       int s=0, c=cQ; ipoint p = P[0].ip;
       while(c--)if(Q[c].rx.mn < p.x && p.x < Q[c].rx.mx)
         {int sgn = 0 < area(p, Q[c].ip, Q[c+1].ip);
           s += sgn != Q[c].ip.x < Q[c+1].ip.x ? 0 : (sgn?-1:1); }
        {int j; for(j=0; j<cP; ++j)
                {if(s) cntrib(P[j].ip, P[j+1].ip, s);
             s += P[j].in;}}}
    inness(ipa, na, ipb, nb); inness(ipb, nb, ipa, na);}
  return s/ascale;}
}

// int main(){point
 /*
   a[] = {{2,3}, {2,3}, {2,3}, {2,4}, {3,3}, {2,3}, {2,3}},
   b[] = {{1,1}, {1,4}, {4,4}, {4,1}, {1,1}}; // 1/2, 1/2
// The redundant vertices above are to provoke errors as good test cases should.
// It is not necessary to duplicate the first vertex at the end.
   a[] = {{1,7}, {4,7}, {4, 6}, {2,6}, {2, 3}, {4,3}, {4,2}, {1,2}},
   b[] = {{3,1}, {5,1}, {5,4}, {3,4}, {3,5}, {6,5}, {6,0}, {3,0}}; // 0, 9
   a[] = {{1,1}, {1,2}, {2,1}, {2,2}},
   b[] = {{0,0}, {0,4}, {4,4}, {4,0}}; // 0, 1/2
   a[] = {{0,0}, {3,0}, {3,2}, {1,2}, {1,1}, {2,1}, {2,3}, {0,3}},
   b[] = {{0,0}, {0,4}, {4,4}, {4,0}}; // -9, 11
   a[] = {{0,0}, {1,0}, {0,1}},
   b[] = {{0,0}, {0,1}, {1,1}, {1,0}}; // -1/2, 1/2
*/
//   a[] = {{1,3}, {2,3}, {2,0}, {1,0}},
//   b[] = {{0,1}, {3,1}, {3,2}, {0,2}}; // -1, 3
/*
   a[] ={{0,0}, {0,2}, {2,2}, {2,0}},
   b[] = {{1, 1}, {3, 1}, {3, 3}, {1, 3}}; // -1, 4
   a[] = {{0,0}, {0,4}, {4,4}, {4,0}},
   b[] = {{1,1}, {1,2}, {2,2}, {2,1}}; // 1, 16
 */
//   {const int ac = sizeof(a)/sizeof(a[0]), bc = sizeof(b)/sizeof(b[0]);
//   printf("%e, %e\n", inter(a, ac, b, bc), inter(a, ac, a, ac)); return 0;}}
