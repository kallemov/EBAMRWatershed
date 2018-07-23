#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <math.h>
#include "RealVect.H"
#include "ParserFunc.H"

#include "NamespaceHeader.H"


ParserFunc::ParserFunc()
    :m_isDefined(false),
     m_muParser(NULL),
     m_point(RealVect(1.2345e+123,1.2345e+123,1.2345e+123)),
     m_depth(1.2345e+123),
     m_time(1.2345e+123),
     m_g(1.2345e+123),
     m_psi(1.2345e+123),
     m_slope(1.2345e+123),
     m_rho(1.2345e+123),
     m_phi(1.2345e+123),
     m_slope_x(1.2345e+123),
     m_slope_y(1.2345e+123),
     m_manncoef(1.2345e+123)
{
}

ParserFunc::ParserFunc(const std::string a_expression)
    :m_isDefined(false),
     m_muParser(NULL),
     m_point(RealVect(1.2345e+123,1.2345e+123,1.2345e+123)),
     m_depth(1.2345e+123),
     m_time(1.2345e+123),
     m_g(1.2345e+123),
     m_psi(1.2345e+123),
     m_slope(1.2345e+123),
     m_rho(1.2345e+123),
     m_phi(1.2345e+123),
     m_slope_x(1.2345e+123),
     m_slope_y(1.2345e+123),
     m_manncoef(1.2345e+123)
{
    define(a_expression);
}

ParserFunc::~ParserFunc()
{
  delete m_muParser;
}

void 
ParserFunc::define(const std::string a_expression)
{
  CH_assert(a_expression != "");

  m_muParser = new mu::Parser();
  m_muParser->SetExpr(a_expression);

  defineConsts();
  defineVars();
  m_isDefined = true;
}

void 
ParserFunc::defineConsts()
{
  m_muParser->DefineConst("pi", M_PI);
  m_muParser->DefineConst("Pi", M_PI);
  m_muParser->DefineConst("PI", M_PI);
}

void 
ParserFunc::defineVars()
{
  m_muParser->DefineVar("T", &m_time);
  m_muParser->DefineVar("t", &m_time);

  m_muParser->DefineVar("X", &m_point[0]);
  m_muParser->DefineVar("x", &m_point[0]);
  m_muParser->DefineVar("Y", &m_point[1]);
  m_muParser->DefineVar("y", &m_point[1]);
  m_muParser->DefineVar("Z", &m_point[2]);
  m_muParser->DefineVar("z", &m_point[2]);
 
  m_muParser->DefineVar("g", &m_g);
  m_muParser->DefineVar("G", &m_g);

  m_muParser->DefineVar("depth", &m_depth);
  m_muParser->DefineVar("psi", &m_psi);
  m_muParser->DefineVar("slope", &m_slope);
  m_muParser->DefineVar("rho", &m_rho); 
  m_muParser->DefineVar("phi", &m_phi);
  m_muParser->DefineVar("slope_x", &m_slope_x);
  m_muParser->DefineVar("slope_y", &m_slope_y);
  m_muParser->DefineVar("manncoef", &m_manncoef);
}

void 
ParserFunc::setPoint(const RealVect& a_point)
{
  m_point = a_point;
}


void 
ParserFunc::setTime(const Real a_time)
{
  m_time = a_time;
}

void 
ParserFunc::setGravitationalConstant(const Real a_g)
{
  m_g = a_g;
}

void 
ParserFunc::setDepth(const Real a_depth)
{
  m_depth = a_depth;
}

void 
ParserFunc::setPressure(const Real a_psi)
{
  m_psi = a_psi;
}

void 
ParserFunc::setSlope(const Real a_slope)
{
  m_slope = a_slope;
}

Real ParserFunc::Eval() const
{
  CH_assert(m_isDefined);
  return m_muParser->Eval();
}

Real 
ParserFunc::Eval(const ParserInput& a_input)
{
  CH_assert(m_isDefined);
  m_point    = a_input.point;
  m_time     = a_input.time;
  m_depth    = a_input.depth;
  m_psi      = a_input.psi;
  m_slope    = a_input.slope;
  m_rho      = a_input.rho; 
  m_phi      = a_input.phi; 
  m_slope_x  = a_input.slope_x;
  m_slope_y  = a_input.slope_y;
  m_manncoef = a_input.manncoef;
    
  return m_muParser->Eval();
}
#include "NamespaceFooter.H"
