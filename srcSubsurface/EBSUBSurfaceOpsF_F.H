#ifndef _EBSUBSURFACEOPSF_F_H_
#define _EBSUBSURFACEOPSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETSATURATIONPOINT 
#define GUARDGETSATURATIONPOINT 
// Prototype for Fortran procedure getsaturationpoint ...
//
void FORTRAN_NAME( GETSATURATIONPOINT ,getsaturationpoint )(
      CHFp_REAL(res)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_REAL(ssat)
      ,CHFp_CONST_REAL(sres)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil) );

#define FORT_GETSATURATIONPOINT FORTRAN_NAME( inlineGETSATURATIONPOINT, inlineGETSATURATIONPOINT)
#define FORTNT_GETSATURATIONPOINT FORTRAN_NAME( GETSATURATIONPOINT, getsaturationpoint)

inline void FORTRAN_NAME(inlineGETSATURATIONPOINT, inlineGETSATURATIONPOINT)(
      CHFp_REAL(res)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_REAL(ssat)
      ,CHFp_CONST_REAL(sres)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil) )
{
 CH_TIMELEAF("FORT_GETSATURATIONPOINT");
 FORTRAN_NAME( GETSATURATIONPOINT ,getsaturationpoint )(
      CHFt_REAL(res)
      ,CHFt_CONST_REAL(h)
      ,CHFt_CONST_REAL(ssat)
      ,CHFt_CONST_REAL(sres)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(nsoil) );
}
#endif  // GUARDGETSATURATIONPOINT 

#ifndef GUARDGETSATURATIONTERMSPOINT 
#define GUARDGETSATURATIONTERMSPOINT 
// Prototype for Fortran procedure getsaturationtermspoint ...
//
void FORTRAN_NAME( GETSATURATIONTERMSPOINT ,getsaturationtermspoint )(
      CHFp_REAL(res)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_REAL(hold)
      ,CHFp_CONST_REAL(ss)
      ,CHFp_CONST_REAL(phi)
      ,CHFp_CONST_REAL(rho)
      ,CHFp_CONST_REAL(ssat)
      ,CHFp_CONST_REAL(sres)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil) );

#define FORT_GETSATURATIONTERMSPOINT FORTRAN_NAME( inlineGETSATURATIONTERMSPOINT, inlineGETSATURATIONTERMSPOINT)
#define FORTNT_GETSATURATIONTERMSPOINT FORTRAN_NAME( GETSATURATIONTERMSPOINT, getsaturationtermspoint)

inline void FORTRAN_NAME(inlineGETSATURATIONTERMSPOINT, inlineGETSATURATIONTERMSPOINT)(
      CHFp_REAL(res)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_REAL(hold)
      ,CHFp_CONST_REAL(ss)
      ,CHFp_CONST_REAL(phi)
      ,CHFp_CONST_REAL(rho)
      ,CHFp_CONST_REAL(ssat)
      ,CHFp_CONST_REAL(sres)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil) )
{
 CH_TIMELEAF("FORT_GETSATURATIONTERMSPOINT");
 FORTRAN_NAME( GETSATURATIONTERMSPOINT ,getsaturationtermspoint )(
      CHFt_REAL(res)
      ,CHFt_CONST_REAL(h)
      ,CHFt_CONST_REAL(hold)
      ,CHFt_CONST_REAL(ss)
      ,CHFt_CONST_REAL(phi)
      ,CHFt_CONST_REAL(rho)
      ,CHFt_CONST_REAL(ssat)
      ,CHFt_CONST_REAL(sres)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(nsoil) );
}
#endif  // GUARDGETSATURATIONTERMSPOINT 

#ifndef GUARDADDSATURATIONTERMS 
#define GUARDADDSATURATIONTERMS 
// Prototype for Fortran procedure addsaturationterms ...
//
void FORTRAN_NAME( ADDSATURATIONTERMS ,addsaturationterms )(
      CHFp_FRA1(fvec)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_FRA1(hold)
      ,CHFp_CONST_FRA1(ss)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(rho)
      ,CHFp_CONST_REAL(ssat)
      ,CHFp_CONST_REAL(sres)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil)
      ,CHFp_BOX(box) );

#define FORT_ADDSATURATIONTERMS FORTRAN_NAME( inlineADDSATURATIONTERMS, inlineADDSATURATIONTERMS)
#define FORTNT_ADDSATURATIONTERMS FORTRAN_NAME( ADDSATURATIONTERMS, addsaturationterms)

inline void FORTRAN_NAME(inlineADDSATURATIONTERMS, inlineADDSATURATIONTERMS)(
      CHFp_FRA1(fvec)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_FRA1(hold)
      ,CHFp_CONST_FRA1(ss)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(rho)
      ,CHFp_CONST_REAL(ssat)
      ,CHFp_CONST_REAL(sres)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADDSATURATIONTERMS");
 FORTRAN_NAME( ADDSATURATIONTERMS ,addsaturationterms )(
      CHFt_FRA1(fvec)
      ,CHFt_CONST_FRA1(h)
      ,CHFt_CONST_FRA1(hold)
      ,CHFt_CONST_FRA1(ss)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(rho)
      ,CHFt_CONST_REAL(ssat)
      ,CHFt_CONST_REAL(sres)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(nsoil)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADDSATURATIONTERMS 

#ifndef GUARDCCPHARMMEANCELLTOFACE 
#define GUARDCCPHARMMEANCELLTOFACE 
// Prototype for Fortran procedure ccpharmmeancelltoface ...
//
void FORTRAN_NAME( CCPHARMMEANCELLTOFACE ,ccpharmmeancelltoface )(
      CHFp_FRA1(faceval)
      ,CHFp_CONST_FRA1(cellval)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) );

#define FORT_CCPHARMMEANCELLTOFACE FORTRAN_NAME( inlineCCPHARMMEANCELLTOFACE, inlineCCPHARMMEANCELLTOFACE)
#define FORTNT_CCPHARMMEANCELLTOFACE FORTRAN_NAME( CCPHARMMEANCELLTOFACE, ccpharmmeancelltoface)

inline void FORTRAN_NAME(inlineCCPHARMMEANCELLTOFACE, inlineCCPHARMMEANCELLTOFACE)(
      CHFp_FRA1(faceval)
      ,CHFp_CONST_FRA1(cellval)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_CCPHARMMEANCELLTOFACE");
 FORTRAN_NAME( CCPHARMMEANCELLTOFACE ,ccpharmmeancelltoface )(
      CHFt_FRA1(faceval)
      ,CHFt_CONST_FRA1(cellval)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDCCPHARMMEANCELLTOFACE 

#ifndef GUARDRELATIVEPERMEABILITY 
#define GUARDRELATIVEPERMEABILITY 
// Prototype for Fortran procedure relativepermeability ...
//
void FORTRAN_NAME( RELATIVEPERMEABILITY ,relativepermeability )(
      CHFp_REAL(kr)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil)
      ,CHFp_CONST_REAL(lambda) );

#define FORT_RELATIVEPERMEABILITY FORTRAN_NAME( inlineRELATIVEPERMEABILITY, inlineRELATIVEPERMEABILITY)
#define FORTNT_RELATIVEPERMEABILITY FORTRAN_NAME( RELATIVEPERMEABILITY, relativepermeability)

inline void FORTRAN_NAME(inlineRELATIVEPERMEABILITY, inlineRELATIVEPERMEABILITY)(
      CHFp_REAL(kr)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil)
      ,CHFp_CONST_REAL(lambda) )
{
 CH_TIMELEAF("FORT_RELATIVEPERMEABILITY");
 FORTRAN_NAME( RELATIVEPERMEABILITY ,relativepermeability )(
      CHFt_REAL(kr)
      ,CHFt_CONST_REAL(h)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(nsoil)
      ,CHFt_CONST_REAL(lambda) );
}
#endif  // GUARDRELATIVEPERMEABILITY 

#ifndef GUARDMULTBYHARMMEANKR 
#define GUARDMULTBYHARMMEANKR 
// Prototype for Fortran procedure multbyharmmeankr ...
//
void FORTRAN_NAME( MULTBYHARMMEANKR ,multbyharmmeankr )(
      CHFp_FRA1(faceval)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil)
      ,CHFp_CONST_REAL(lambda)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) );

#define FORT_MULTBYHARMMEANKR FORTRAN_NAME( inlineMULTBYHARMMEANKR, inlineMULTBYHARMMEANKR)
#define FORTNT_MULTBYHARMMEANKR FORTRAN_NAME( MULTBYHARMMEANKR, multbyharmmeankr)

inline void FORTRAN_NAME(inlineMULTBYHARMMEANKR, inlineMULTBYHARMMEANKR)(
      CHFp_FRA1(faceval)
      ,CHFp_CONST_FRA1(h)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(nsoil)
      ,CHFp_CONST_REAL(lambda)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_MULTBYHARMMEANKR");
 FORTRAN_NAME( MULTBYHARMMEANKR ,multbyharmmeankr )(
      CHFt_FRA1(faceval)
      ,CHFt_CONST_FRA1(h)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(nsoil)
      ,CHFt_CONST_REAL(lambda)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDMULTBYHARMMEANKR 

#ifndef GUARDMACDIVERGEF 
#define GUARDMACDIVERGEF 
// Prototype for Fortran procedure macdivergef ...
//
void FORTRAN_NAME( MACDIVERGEF ,macdivergef )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) );

#define FORT_MACDIVERGEF FORTRAN_NAME( inlineMACDIVERGEF, inlineMACDIVERGEF)
#define FORTNT_MACDIVERGEF FORTRAN_NAME( MACDIVERGEF, macdivergef)

inline void FORTRAN_NAME(inlineMACDIVERGEF, inlineMACDIVERGEF)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_MACDIVERGEF");
 FORTRAN_NAME( MACDIVERGEF ,macdivergef )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(divf)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDMACDIVERGEF 

#ifndef GUARDCCPAVEFACETOCELL 
#define GUARDCCPAVEFACETOCELL 
// Prototype for Fortran procedure ccpavefacetocell ...
//
void FORTRAN_NAME( CCPAVEFACETOCELL ,ccpavefacetocell )(
      CHFp_FRA1(cellgrad)
      ,CHFp_CONST_FRA1(facegrad)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(cellbox) );

#define FORT_CCPAVEFACETOCELL FORTRAN_NAME( inlineCCPAVEFACETOCELL, inlineCCPAVEFACETOCELL)
#define FORTNT_CCPAVEFACETOCELL FORTRAN_NAME( CCPAVEFACETOCELL, ccpavefacetocell)

inline void FORTRAN_NAME(inlineCCPAVEFACETOCELL, inlineCCPAVEFACETOCELL)(
      CHFp_FRA1(cellgrad)
      ,CHFp_CONST_FRA1(facegrad)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_CCPAVEFACETOCELL");
 FORTRAN_NAME( CCPAVEFACETOCELL ,ccpavefacetocell )(
      CHFt_FRA1(cellgrad)
      ,CHFt_CONST_FRA1(facegrad)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDCCPAVEFACETOCELL 

}

#endif
