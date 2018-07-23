#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RealVect.H"
#include "ParserBCValue.H"

#include "NamespaceHeader.H"
ParserBCValue::ParserBCValue()
{
  //m_gconst=1.0;
  m_isDefined = false;
}

ParserBCValue::~ParserBCValue()
{
}

void ParserBCValue::define(const RefCountedPtr<ParserFunc> a_muParserFuncPtr, const int a_BCtype)
{
  m_isDefined = true;
  m_muParserFuncPtr = a_muParserFuncPtr;
  m_BCtype = a_BCtype;
}

Real ParserBCValue::value(const RealVect& a_point,
			   const RealVect& a_normal,
			   const Real&     a_time,
			   const int&      a_comp) const
{
  CH_assert(m_isDefined);

  m_muParserFuncPtr->setPoint(a_point);
  m_muParserFuncPtr->setTime(a_time);

  Real value = m_muParserFuncPtr->Eval();
  
  //add hydsrostatic term
  if (m_BCtype==1)
    {
      //      value += a_point[2]*m_gconst;
      value += a_point[2];
    } 
  //change the sign because it changes in divergence operator
  else if (m_BCtype==0)
  {
    value *=-1.0;
  }
 
  return value;
}

Real ParserBCValue::parserValue(const RealVect& a_point,
				const RealVect& a_normal,
				const Real&     a_time,
				const int&      a_comp) const
{
  CH_assert(m_isDefined);

  m_muParserFuncPtr->setPoint(a_point);
  m_muParserFuncPtr->setTime(a_time);
  Real value = m_muParserFuncPtr->Eval();
  
  return value;
}

// void ParserBCValue::setGravityConstant(Real a_gconst)
// {
//   m_gconst=a_gconst;
//   m_muParserFuncPtr->setGravitationalConstant(m_gconst);
// }
#include "NamespaceFooter.H"
