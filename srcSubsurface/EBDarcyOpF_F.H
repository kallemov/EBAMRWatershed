#ifndef _EBDARCYOPF_F_H_
#define _EBDARCYOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDDARCYEBCOREGAPPLYDOMAINFLUX 
#define GUARDDARCYEBCOREGAPPLYDOMAINFLUX 
// Prototype for Fortran procedure darcyebcoregapplydomainflux ...
//
void FORTRAN_NAME( DARCYEBCOREGAPPLYDOMAINFLUX ,darcyebcoregapplydomainflux )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_FRA1(bc)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_DARCYEBCOREGAPPLYDOMAINFLUX FORTRAN_NAME( inlineDARCYEBCOREGAPPLYDOMAINFLUX, inlineDARCYEBCOREGAPPLYDOMAINFLUX)
#define FORTNT_DARCYEBCOREGAPPLYDOMAINFLUX FORTRAN_NAME( DARCYEBCOREGAPPLYDOMAINFLUX, darcyebcoregapplydomainflux)

inline void FORTRAN_NAME(inlineDARCYEBCOREGAPPLYDOMAINFLUX, inlineDARCYEBCOREGAPPLYDOMAINFLUX)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_FRA1(bc)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DARCYEBCOREGAPPLYDOMAINFLUX");
 FORTRAN_NAME( DARCYEBCOREGAPPLYDOMAINFLUX ,darcyebcoregapplydomainflux )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(faceflux)
      ,CHFt_CONST_FRA1(bc)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDARCYEBCOREGAPPLYDOMAINFLUX 

#ifndef GUARDDARCYCONDUCTIVITYINPLACE 
#define GUARDDARCYCONDUCTIVITYINPLACE 
// Prototype for Fortran procedure darcyconductivityinplace ...
//
void FORTRAN_NAME( DARCYCONDUCTIVITYINPLACE ,darcyconductivityinplace )(
      CHFp_FRA1(opphidir)
      ,CHFp_FRA1(darcyflux0)
      ,CHFp_FRA1(darcyflux1)
      ,CHFp_FRA1(darcyflux2)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) );

#define FORT_DARCYCONDUCTIVITYINPLACE FORTRAN_NAME( inlineDARCYCONDUCTIVITYINPLACE, inlineDARCYCONDUCTIVITYINPLACE)
#define FORTNT_DARCYCONDUCTIVITYINPLACE FORTRAN_NAME( DARCYCONDUCTIVITYINPLACE, darcyconductivityinplace)

inline void FORTRAN_NAME(inlineDARCYCONDUCTIVITYINPLACE, inlineDARCYCONDUCTIVITYINPLACE)(
      CHFp_FRA1(opphidir)
      ,CHFp_FRA1(darcyflux0)
      ,CHFp_FRA1(darcyflux1)
      ,CHFp_FRA1(darcyflux2)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DARCYCONDUCTIVITYINPLACE");
 FORTRAN_NAME( DARCYCONDUCTIVITYINPLACE ,darcyconductivityinplace )(
      CHFt_FRA1(opphidir)
      ,CHFt_FRA1(darcyflux0)
      ,CHFt_FRA1(darcyflux1)
      ,CHFt_FRA1(darcyflux2)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(b0)
      ,CHFt_CONST_FRA1(b1)
      ,CHFt_CONST_FRA1(b2)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDARCYCONDUCTIVITYINPLACE 

#ifndef GUARDDARCYINCRAPPLYEBCO 
#define GUARDDARCYINCRAPPLYEBCO 
// Prototype for Fortran procedure darcyincrapplyebco ...
//
void FORTRAN_NAME( DARCYINCRAPPLYEBCO ,darcyincrapplyebco )(
      CHFp_FRA1(lhs)
      ,CHFp_CONST_FRA1(interiorflux)
      ,CHFp_CONST_FRA1(domainfluxlo)
      ,CHFp_CONST_FRA1(domainfluxhi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(facedir) );

#define FORT_DARCYINCRAPPLYEBCO FORTRAN_NAME( inlineDARCYINCRAPPLYEBCO, inlineDARCYINCRAPPLYEBCO)
#define FORTNT_DARCYINCRAPPLYEBCO FORTRAN_NAME( DARCYINCRAPPLYEBCO, darcyincrapplyebco)

inline void FORTRAN_NAME(inlineDARCYINCRAPPLYEBCO, inlineDARCYINCRAPPLYEBCO)(
      CHFp_FRA1(lhs)
      ,CHFp_CONST_FRA1(interiorflux)
      ,CHFp_CONST_FRA1(domainfluxlo)
      ,CHFp_CONST_FRA1(domainfluxhi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_DARCYINCRAPPLYEBCO");
 FORTRAN_NAME( DARCYINCRAPPLYEBCO ,darcyincrapplyebco )(
      CHFt_FRA1(lhs)
      ,CHFt_CONST_FRA1(interiorflux)
      ,CHFt_CONST_FRA1(domainfluxlo)
      ,CHFt_CONST_FRA1(domainfluxhi)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(lobox)
      ,CHFt_BOX(hibox)
      ,CHFt_BOX(centerbox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDDARCYINCRAPPLYEBCO 

#ifndef GUARDDARCYDECRINVRELCOEFEBCO 
#define GUARDDARCYDECRINVRELCOEFEBCO 
// Prototype for Fortran procedure darcydecrinvrelcoefebco ...
//
void FORTRAN_NAME( DARCYDECRINVRELCOEFEBCO ,darcydecrinvrelcoefebco )(
      CHFp_FRA1(relcoef)
      ,CHFp_FRA1(bcoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir) );

#define FORT_DARCYDECRINVRELCOEFEBCO FORTRAN_NAME( inlineDARCYDECRINVRELCOEFEBCO, inlineDARCYDECRINVRELCOEFEBCO)
#define FORTNT_DARCYDECRINVRELCOEFEBCO FORTRAN_NAME( DARCYDECRINVRELCOEFEBCO, darcydecrinvrelcoefebco)

inline void FORTRAN_NAME(inlineDARCYDECRINVRELCOEFEBCO, inlineDARCYDECRINVRELCOEFEBCO)(
      CHFp_FRA1(relcoef)
      ,CHFp_FRA1(bcoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir) )
{
 CH_TIMELEAF("FORT_DARCYDECRINVRELCOEFEBCO");
 FORTRAN_NAME( DARCYDECRINVRELCOEFEBCO ,darcydecrinvrelcoefebco )(
      CHFt_FRA1(relcoef)
      ,CHFt_FRA1(bcoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box)
      ,CHFt_REAL(dx)
      ,CHFt_INT(idir) );
}
#endif  // GUARDDARCYDECRINVRELCOEFEBCO 

#ifndef GUARDDARCYINVERTLAMBDAEBCO 
#define GUARDDARCYINVERTLAMBDAEBCO 
// Prototype for Fortran procedure darcyinvertlambdaebco ...
//
void FORTRAN_NAME( DARCYINVERTLAMBDAEBCO ,darcyinvertlambdaebco )(
      CHFp_FRA1(lambda)
      ,CHFp_REAL(safety)
      ,CHFp_BOX(box) );

#define FORT_DARCYINVERTLAMBDAEBCO FORTRAN_NAME( inlineDARCYINVERTLAMBDAEBCO, inlineDARCYINVERTLAMBDAEBCO)
#define FORTNT_DARCYINVERTLAMBDAEBCO FORTRAN_NAME( DARCYINVERTLAMBDAEBCO, darcyinvertlambdaebco)

inline void FORTRAN_NAME(inlineDARCYINVERTLAMBDAEBCO, inlineDARCYINVERTLAMBDAEBCO)(
      CHFp_FRA1(lambda)
      ,CHFp_REAL(safety)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DARCYINVERTLAMBDAEBCO");
 FORTRAN_NAME( DARCYINVERTLAMBDAEBCO ,darcyinvertlambdaebco )(
      CHFt_FRA1(lambda)
      ,CHFt_REAL(safety)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDARCYINVERTLAMBDAEBCO 

#ifndef GUARDDARCYGETFLUXEBCO 
#define GUARDDARCYGETFLUXEBCO 
// Prototype for Fortran procedure darcygetfluxebco ...
//
void FORTRAN_NAME( DARCYGETFLUXEBCO ,darcygetfluxebco )(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(opphibox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir) );

#define FORT_DARCYGETFLUXEBCO FORTRAN_NAME( inlineDARCYGETFLUXEBCO, inlineDARCYGETFLUXEBCO)
#define FORTNT_DARCYGETFLUXEBCO FORTRAN_NAME( DARCYGETFLUXEBCO, darcygetfluxebco)

inline void FORTRAN_NAME(inlineDARCYGETFLUXEBCO, inlineDARCYGETFLUXEBCO)(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(opphibox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_DARCYGETFLUXEBCO");
 FORTRAN_NAME( DARCYGETFLUXEBCO ,darcygetfluxebco )(
      CHFt_FRA1(flux)
      ,CHFt_CONST_FRA1(bcoef)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_BOX(opphibox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDDARCYGETFLUXEBCO 

}

#endif
