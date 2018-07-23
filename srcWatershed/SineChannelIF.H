#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SINECHANNELIF_H_
#define _SINECHANNELIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a 
    sine channel
    f(x,y,z) = z - z0 + depth*sin(pi/2+pi*(x-x0)/width)+ slope*(y-y0)
    
 */
class SineChannelIF: public BaseIF
{
public:
  SineChannelIF(const RealVect & a_point,
		const Real     & a_width,
		const Real     & a_slope,
		const Real     & a_depth,
		const bool       a_inside = false);
    
  /// Destructor
  virtual ~SineChannelIF()
  {
  }


  virtual Real value(const RealVect& a_point) const;
    
  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
		     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:

  RealVect m_point;
  Real     m_width;      
  Real     m_slope;      
  Real     m_depth;      
  bool     m_inside;
private:
  SineChannelIF();

};

#include "NamespaceFooter.H"
#endif
