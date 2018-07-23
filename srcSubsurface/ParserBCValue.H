#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PARSERBCVALUE_H_
#define _PARSERBCVALUE_H_

#include "RealVect.H"

#include "BaseBCValue.H"
#include "ParserFunc.H"
#include "MayDay.H"
#include "NamespaceHeader.H"

class ParserBCValue: public BaseBCValue
{
public:
  ParserBCValue();

  virtual ~ParserBCValue();

  virtual void define(const RefCountedPtr<ParserFunc> a_muParserFuncPtr, const int a_BCtype);

  virtual Real value(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;

  Real parserValue(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;
  //void setGravityConstant(Real a_gconst);
private:
  bool m_isDefined;
  int m_BCtype;
  RefCountedPtr<ParserFunc> m_muParserFuncPtr;
  //Real m_gconst;
};

#include "NamespaceFooter.H"
#endif
