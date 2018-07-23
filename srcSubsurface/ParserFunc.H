#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PARSERFUNC_H_
#define _PARSERFUNC_H_

#include "muParser.h"
#include "RealVect.H"
#include "MayDay.H"

#include "NamespaceHeader.H"

typedef struct ParserInput ParserInput;
struct ParserInput{
    RealVect point;
    Real     time;
    Real     depth;
    Real     psi;
    Real     slope;
    Real     rho;
    Real     phi;
    Real     slope_x;
    Real     slope_y;
    Real     manncoef;
    
    ParserInput():
	point(RealVect(1.2345e+123,1.2345e+123,1.2345e+123)),
	time(1.2345e+123),
	depth(1.2345e+123),
	psi(1.2345e+123),
	slope(1.2345e+123),
        rho(1.2345e+123),
	phi(1.2345e+123),
	slope_x(1.2345e+123),
	slope_y(1.2345e+123),
	manncoef(1.2345e+123){};
};

class ParserFunc
{
public:
    ParserFunc();
    ParserFunc(const std::string a_expression);
    
    virtual ~ParserFunc();
    
    virtual void define(const std::string a_expression);
    virtual void setPoint(const RealVect& a_point);
    virtual void setTime(const Real a_time);
    virtual void setGravitationalConstant(const Real a_g);
    virtual void setDepth(const Real a_depth);
    virtual void setPressure(const Real a_psi);
    virtual void setSlope(const Real a_slope);
    
    // virtual void setDensity(const Real a_rho);
    // virtual void setVelocity(const RealVect& a_vel);
    // virtual void setPressure(const RealVect& a_vel);
    
    virtual Real Eval() const;
    virtual Real Eval(const ParserInput& a_input);
    
protected:
    
    virtual void defineConsts();
    virtual void defineVars();
    
private:
    bool m_isDefined;
    mu::Parser* m_muParser;
    
    RealVect m_point; //coordinate
    Real     m_depth; //depth from the surface
    Real     m_time; // time
    Real     m_g;   //gravitational constant
    Real     m_psi; //preasurehead
    Real     m_slope; // slope for current direction
    Real     m_rho; //density
    Real     m_phi; //porosity
    Real     m_slope_x;//slope in x direction
    Real     m_slope_y;//slope in y direction
    Real     m_manncoef; //Manning coefficient
};
#include "NamespaceFooter.H"
#endif
