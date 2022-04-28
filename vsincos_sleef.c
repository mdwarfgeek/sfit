/* Implementation of v_sincos using SIMD functions from the SLEEF library. */

void FUNC (double v, double *t, int ndp,
           double *sinarg, double *cosarg) {
  int idp, iex;
  double x;

  REGD v_v, v_t, v_x, v_i;
  REGD2 v_sc;
  Sleef_double2 sc;

  /* How much can we compute using packed doubles? */
  iex = VLEN * (ndp / VLEN);

  /* Packed double calculation */
  v_v = MM_PREFIX(set1_pd)(v);
  
  for(idp = 0; idp < iex; idp += VLEN) {
    /* Load elements from t */
    v_t = MM_PREFIX(loadu_pd)(&(t[idp]));
    
    /* x = v * t */
    v_x = MM_PREFIX(mul_pd)(v_v, v_t);
    
    /* x -= rint(x) */
    /* ROUNDPD is an SSE4.1 instruction, and this line is where the
       dependency comes from.  It could potentially be avoided by
       doing the "convert to integer and back again" trick but the
       range of arguments this can cope with is limited to the range
       of a 32-bit integer, which doesn't really solve the problem
       we're trying to solve by range reducing (see below).  So for
       now we require SSE4.1 to be able to use ROUNDPD to do it
       properly. */
    v_i = MM_PREFIX(round_pd)(v_x, 4);
    v_x = MM_PREFIX(sub_pd)(v_x, v_i);
    
    /* 2*x */
    v_x = MM_PREFIX(add_pd)(v_x, v_x);

    /* Call sincospi.  Allowed range of arguments is [-1e9, 1e9]
       so there's an explicit range reduction step above to make
       sure this is always true. */
    v_sc = SLEEF_FUNC(Sleef_sincospi, u35)(v_x);

    /* Store results */
    MM_PREFIX(storeu_pd)(&(sinarg[idp]), v_sc.x);
    MM_PREFIX(storeu_pd)(&(cosarg[idp]), v_sc.y);
  }

  /* Do the rest as scalars in plain C */
  for(idp = iex; idp < ndp; idp++) {
    /* Range reduction */
    x = v * t[idp];
    x -= rint(x);

    /* 2*x */
    x += x;

    /* Call sincospi */
    sc = Sleef_sincospi_u35(x);

    /* Store results */
    sinarg[idp] = sc.x;
    cosarg[idp] = sc.y;
  }
}
