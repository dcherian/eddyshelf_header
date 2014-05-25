/* Solver */
#define SOLVE3D

/*#define SPLINES  - apparently problematic */

/* Analytic Equations */

/* Boundary Conditions */
#define ANA_BTFLUX
#define ANA_STFLUX
#define ANA_SMFLUX /* surface momentum flux*/
#define ANA_SSFLUX
#define ANA_BSFLUX

#define UV_ADV
#define UV_COR
#define RADIATION_2D
#define DJ_GRADPS
#define PERFECT_RESTART
/*#define N2S2_HORAVG*/
#undef  NONLIN_EOS
#undef  SALINITY

/* FLOATS */
#define FLOATS
#define FLOAT_VWALK
#define FLOAT_STICKY

#define DC_FLOATS_DEPLOYMENT

/* Advection Schemes */

#undef SPLIT_SCHEME

#ifdef SPLIT_SCHEME

#define UV_U3ADV_SPLIT
#define TS_U3ADV_SPLIT
#define UV_VIS4
#define TS_DIF4

#else

/*#define TS_A4HADVECTION
  #define TS_A4VADVECTION*/
#define TS_MPDATA /* Horizontal advection for T,S */
#define UV_C3HADVECTION

#define VISC_GRID
#define DIFF_GRID

#define UV_VIS2
#define TS_DIF2
#define UV_VIS4
#define TS_DIF4
#endif

/*#define UV_SVADVECTION - for shallow, well resolved domains*/

/* Viscosity & sponges */
#define UV_LDRAG
#define LIMIT_BSTRESS

/*#define ANA_DRAG
#define UV_DRAG_GRID

#define DC_SPRDRG
#define DC_SPRDRG_EAST
#define DC_SPRDRG_WEST
#define DC_SPRDRG_NORTH*/

#define ANA_SPONGE
#define DC_SPONGE
#define DC_SPONGE_WEST
#define DC_SPONGE_EAST
#define DC_SPONGE_NORTH
/*#define DC_SPONGE_SOUTH*/

/*#define TS_FIXED - diagnostic spinup*/

/* Nudging options */
/*#define M2CLIMATOLOGY
#define M2CLM_NUDGING

#define DC_NUDGING
#define DC_NUDGE_NORTH
#define DC_NUDGE_WEST
#define DC_NUDGE_EAST
/*#define DC_NUDGE_SOUTH*/

/* Passive Tracer */
#define T_PASSIVE
#define ANA_SPFLUX
#define ANA_BPFLUX 

/* Mixing */
#define MIX_GEO_UV4
#define MIX_GEO_TS4

#define MIX_S_UV2
#define MIX_GEO_TS2
/*#define MY25_MIXING
  #define KANTHA_CLAYSON*/
#define GLS_MIXING
#define CANUTO_A

/* Outputs */
#define AVERAGES
/*#define DIAGNOSTICS_UV
#define DIAGNOSTICS_TS*/

/* NETCDF 4 */
#define HDF5
#define DEFLATE
