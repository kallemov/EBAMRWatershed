#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SineChannelIF.H"
#include "CONSTANTS.H"
#include "NamespaceHeader.H"

SineChannelIF::
SineChannelIF(const RealVect & a_point,
	      const Real     & a_width,
	      const Real     & a_slope,
	      const Real     & a_depth,
	      const bool       a_inside)
    :BaseIF()
{
    m_point    = a_point;
    m_width    = a_width;
    m_slope    = a_slope;
    m_depth    = a_depth;
    m_inside   = a_inside;
}

Real 
SineChannelIF::
value(const RealVect& a_point) const
{
    
    Real retval = -(a_point[2] - m_point[2]) + m_slope*(a_point[0] - m_point[0]);

    Real yarg = 0.5 + (a_point[1] - m_point[1])/m_width;
    //if (yarg<1.0 && yarg>0.0)
    retval -=  m_depth*sin(yarg*PI);

    if(!m_inside)
	retval = -retval;
    
    return retval;
}

Real
SineChannelIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
		     const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Real retval= LARGEREALVAL;
  int derivativeOrder = a_partialDerivative.sum();
  
  if (derivativeOrder == 0)
  {
      retval = value(a_point);
  }
  else if (derivativeOrder == 1)
  {
      Real yarg = 0.5 + (a_point[1] - m_point[1])/m_width;
      
      retval = (Real)a_partialDerivative[0]*m_slope
	  -(Real)a_partialDerivative[1]*m_depth*PI/m_width*cos(yarg*PI)
	  -(Real)a_partialDerivative[2];
  }
  else
  {
      if (a_partialDerivative[1] == derivativeOrder)
      {
	  Real yarg = 0.5 + (a_point[1] - m_point[1])/m_width;
	  retval = -m_depth*pow(PI/m_width,derivativeOrder)*sin(derivativeOrder*PI/2.0+yarg*PI);
      }
      else
      {
          // mixed partials = 0.0
          retval = 0.0;
      }
  }

  // Change the sign to change inside to outside
  if (!m_inside && derivativeOrder > 0)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* 
SineChannelIF::
newImplicitFunction() const
{
    SineChannelIF* sinePtr = new SineChannelIF(m_point, m_width, m_slope, m_depth,m_inside);
  return static_cast<BaseIF*>(sinePtr);
}

#include "NamespaceFooter.H"
